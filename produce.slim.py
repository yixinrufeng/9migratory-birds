#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate SLiM 5.x nonWF models for bird conservation simulations.

Current version:
- burn-in phase is regulated around ancestral_N
- current phase is regulated around current_N
- future phase starts from the current population after the current phase
- future scenarios set only an upper carrying-capacity ceiling:
  0.25x, 0.5x, 1x, 2x, 5x, and 10x current_N
- future phase has no lower-bound rescue: when N is below the ceiling,
  pop.fitnessScaling remains 1.0
- adult breeding probability is drawn from N(0.8, 0.1) truncated to [0, 1]
- each egg has a 0.39 probability of becoming a recruited offspring
- neutral mutations are simulated
- neutral : deleterious mutation ratio = 1 : 2.31
- deleterious nonsynonymous mutations follow the Kyriazis-style DFE
- individuals are age-structured; one SLiM cycle represents one year
- annual age-specific survival is applied via individual fitnessScaling

Run:
    python3 generate_bird_slim5_models_K_ceiling.py
    bash run_all_slim.sh
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List
import textwrap


@dataclass(frozen=True)
class SpeciesParams:
    species: str
    age_first_repro: int
    current_N: int
    clutch_mean: float
    clutch_sd: float
    max_lifespan: int
    ancestral_N: int
    mutation_rate: float
    generation_time: float


SPECIES_TABLE: List[SpeciesParams] = [
    SpeciesParams(
        species="spoon",
        age_first_repro=3,
        current_N=490,
        clutch_mean=4.0,
        clutch_sd=1.0,
        max_lifespan=22,
        ancestral_N=2982,
        mutation_rate=1.7e-9,
        generation_time=4.5,
    ),
]


N_GENES = 15_000
N_CHROMOSOMES = 44
GENE_LENGTH = 2_000
RECOMBINATION_RATE = 1.0e-9

NEUTRAL_WEIGHT = 1.0
DELETERIOUS_TO_NEUTRAL_RATIO = 2.31

ADULT_REPRODUCTION_PROBABILITY_MEAN = 0.8
ADULT_REPRODUCTION_PROBABILITY_SD = 0.1
CLUTCH_TO_ADULT_PROBABILITY = 0.39
DISASTER_ALPHA = 0.5
DISASTER_BETA = 8.0

CURRENT_PHASE_GENERATIONS = 10
FUTURE_MAX_GENERATIONS = 5_000

ROH_CUTOFF_BP = 1_000_000
LOG_INTERVAL = 100
BURNIN_PROGRESS_INTERVAL = 100
STATS_SAMPLE_SIZE = 100

DFE_MEAN_ABS_S = 0.0131
DFE_SHAPE = 0.186
LETHAL_FRACTION_OF_DELETERIOUS = 0.005

DELETERIOUS_CATEGORY_WEIGHTS = {
    "weak": 1.1273361920,
    "moderate": 0.5693702773,
    "strong": 0.5419366787,
    "very_strong": 0.0598068174,
    "lethal": 0.0115500000,
}

# Multipliers define only the future carrying-capacity ceiling.
# They do not resize the population at future generation 0.
SCENARIOS: Dict[str, float] = {
    "K_0p25x": 0.25,
    "K_0p5x": 0.5,
    "K_1x": 1.0,
    "K_2x": 2.0,
    "K_5x": 5.0,
    "K_10x": 10.0,
    "K_20x": 20.0,
}

MODEL_DIR = Path("slim_models")
RESULT_DIR = Path("slim_results")
RUN_SCRIPT = Path("run_all_slim.sh")


def sanitize_name(name: str) -> str:
    out = []
    for ch in name.strip():
        if ch.isalnum() or ch in {"_", "-"}:
            out.append(ch)
        elif ch.isspace():
            out.append("_")
    return "".join(out) or "species"


def chromosome_lengths() -> List[int]:
    """Distribute all N_GENES genes across N_CHROMOSOMES chromosomes."""
    base = N_GENES // N_CHROMOSOMES
    remainder = N_GENES % N_CHROMOSOMES

    lengths = []
    for i in range(N_CHROMOSOMES):
        genes_on_chr = base + (1 if i < remainder else 0)
        lengths.append(genes_on_chr * GENE_LENGTH)

    assert sum(lengths) == N_GENES * GENE_LENGTH
    return lengths


def life_table(max_lifespan: int) -> List[float]:
    """Age-specific mortality probabilities."""
    if max_lifespan < 1:
        raise ValueError("max_lifespan must be >= 1")

    values = [0.0] * (max_lifespan + 1)
    values[0] = 0.5
    if max_lifespan >= 1:
        values[1] = 0.2

    for age, mortality in [
        (max_lifespan - 4, 0.1),
        (max_lifespan - 3, 0.25),
        (max_lifespan - 2, 0.5),
        (max_lifespan - 1, 0.75),
        (max_lifespan, 1.0),
    ]:
        if 0 <= age <= max_lifespan:
            values[age] = mortality

    return values


def eidos_vector(values: List[float | int]) -> str:
    return "c(" + ", ".join(str(v) for v in values) + ")"


def scaled_integer(value: float) -> int:
    """
    Convert a scaled population size or K ceiling to an integer.
    Uses half-up rounding, so 50 * 0.25 = 12.5 becomes 13, not 12.
    """
    return max(1, int(value + 0.5))


def build_slim_script(
    sp: SpeciesParams,
    scenario_name: str,
    scenario_multiplier: float,
) -> str:
    safe_species = sanitize_name(sp.species)

    future_K_ceiling = scaled_integer(sp.current_N * scenario_multiplier)
    burnin_generations = 10 * sp.ancestral_N

    future_start_tick = burnin_generations + CURRENT_PHASE_GENERATIONS + 1
    end_tick = burnin_generations + CURRENT_PHASE_GENERATIONS + FUTURE_MAX_GENERATIONS

    chrom_lengths = chromosome_lengths()
    genome_length = sum(chrom_lengths)
    life = life_table(sp.max_lifespan)
    out_file = RESULT_DIR / f"{safe_species}.{scenario_name}.tsv"

    template = f'''
// Generated by generate_bird_slim5_models_K_ceiling.py
// Species: {sp.species}
// Scenario: {scenario_name}
// Future K ceiling = {future_K_ceiling}
// Future phase uses an upper carrying-capacity ceiling only.
// There is no lower-bound rescue when N is below the ceiling.

initialize() {{
    initializeSLiMModelType("nonWF");
    initializeSex();

    defineConstant("SPECIES_NAME", "{sp.species}");
    defineConstant("SCENARIO", "{scenario_name}");

    defineConstant("AGE_FIRST_REPRO", {sp.age_first_repro});
    defineConstant("CURRENT_N", {sp.current_N});
    defineConstant("ANCESTRAL_N", {sp.ancestral_N});
    defineConstant("FUTURE_K_CEILING", {future_K_ceiling});

    defineConstant("CLUTCH_MEAN", {sp.clutch_mean});
    defineConstant("CLUTCH_SD", {sp.clutch_sd});
    defineConstant("MAX_LIFESPAN", {sp.max_lifespan});
    defineConstant("GENERATION_TIME", {sp.generation_time});

    defineConstant("BURNIN_GENERATIONS", {burnin_generations});
    defineConstant("CURRENT_PHASE_GENERATIONS", {CURRENT_PHASE_GENERATIONS});
    defineConstant("FUTURE_MAX_GENERATIONS", {FUTURE_MAX_GENERATIONS});
    defineConstant("LOG_INTERVAL", {LOG_INTERVAL});
    defineConstant("BURNIN_PROGRESS_INTERVAL", {BURNIN_PROGRESS_INTERVAL});
    defineConstant("STATS_SAMPLE_SIZE", {STATS_SAMPLE_SIZE});

    defineConstant("N_GENES", {N_GENES});
    defineConstant("N_CHROMOSOMES", {N_CHROMOSOMES});
    defineConstant("GENE_LENGTH", {GENE_LENGTH});
    defineConstant("CHR_LENGTHS", {eidos_vector(chrom_lengths)});
    defineConstant("GENOME_LENGTH", {genome_length});

    defineConstant("MUTATION_RATE", {sp.mutation_rate});
    defineConstant("RECOMBINATION_RATE", {RECOMBINATION_RATE});
    defineConstant("ROH_CUTOFF_BP", {ROH_CUTOFF_BP});

    defineConstant("ADULT_REPRO_PROB_MEAN", {ADULT_REPRODUCTION_PROBABILITY_MEAN});
    defineConstant("ADULT_REPRO_PROB_SD", {ADULT_REPRODUCTION_PROBABILITY_SD});
    defineConstant("CLUTCH_TO_ADULT_PROB", {CLUTCH_TO_ADULT_PROBABILITY});
    defineConstant("DISASTER_ALPHA", {DISASTER_ALPHA});
    defineConstant("DISASTER_BETA", {DISASTER_BETA});
    defineConstant("LIFE_TABLE", {eidos_vector(life)});

    defineConstant("DFE_MEAN_ABS_S", {DFE_MEAN_ABS_S});
    defineConstant("DFE_SHAPE", {DFE_SHAPE});
    defineConstant("OUT_FILE", "{out_file.as_posix()}");

    initializeMutationType("m1", 0.50, "f", 0.0);
    initializeMutationType("m2", 0.45, "f", -0.0001);
    initializeMutationType("m3", 0.20, "f", -0.005);
    initializeMutationType("m4", 0.05, "f", -0.05);
    initializeMutationType("m5", 0.00, "f", -0.5);
    initializeMutationType("m6", 0.00, "f", -1.0);

    m1.convertToSubstitution = T;
    m2.convertToSubstitution = F;
    m3.convertToSubstitution = F;
    m4.convertToSubstitution = F;
    m5.convertToSubstitution = F;
    m6.convertToSubstitution = F;

    initializeGenomicElementType(
        "g1",
        c(m1, m2, m3, m4, m5, m6),
        c({NEUTRAL_WEIGHT},
          {DELETERIOUS_CATEGORY_WEIGHTS['weak']},
          {DELETERIOUS_CATEGORY_WEIGHTS['moderate']},
          {DELETERIOUS_CATEGORY_WEIGHTS['strong']},
          {DELETERIOUS_CATEGORY_WEIGHTS['very_strong']},
          {DELETERIOUS_CATEGORY_WEIGHTS['lethal']})
    );

    for (chrID in 1:N_CHROMOSOMES) {{
        len = CHR_LENGTHS[chrID - 1];
        initializeChromosome(chrID, len, type="A", symbol="chr" + chrID);
        initializeGenomicElement(g1, 0, len - 1);
        initializeMutationRate(MUTATION_RATE);
        initializeRecombinationRate(RECOMBINATION_RATE);
    }}
}}

function (float$)drawTruncatedGammaAbsS(float$ lo, float$ hi)
{{
    while (T) {{
        x = rgamma(1, DFE_MEAN_ABS_S, DFE_SHAPE);
        if ((x >= lo) & (x < hi))
            return x;
    }}

    return NAN;
}}

mutation(m2) {{
    mut.setSelectionCoeff(-drawTruncatedGammaAbsS(0.0, 0.001));
    return T;
}}

mutation(m3) {{
    mut.setSelectionCoeff(-drawTruncatedGammaAbsS(0.001, 0.01));
    return T;
}}

mutation(m4) {{
    mut.setSelectionCoeff(-drawTruncatedGammaAbsS(0.01, 0.1));
    return T;
}}

mutation(m5) {{
    mut.setSelectionCoeff(-drawTruncatedGammaAbsS(0.1, 1.0));
    return T;
}}

function (float)ageSurvivalForIndividuals(object<Individual> inds)
{{
    // Annual age-specific survival probability for each individual.
    // In this age-structured nonWF model, one SLiM cycle represents one year.
    // Offspring are age 0 by default, and individual age increases automatically
    // by one at the end of each cycle.
    if (size(inds) == 0)
        return c();

    ages = inds.age;
    ages[ages > MAX_LIFESPAN] = MAX_LIFESPAN;

    mortality = LIFE_TABLE[ages];
    return 1.0 - mortality;
}}

function (float$)meanAgeSurvivalForIndividuals(object<Individual> inds)
{{
    if (size(inds) == 0)
        return NAN;

    return mean(ageSurvivalForIndividuals(inds));
}}

function (void)applyRegulatedFitnessScaling(object<Subpopulation>$ pop, numeric$ targetN)
{{
    inds = pop.individuals;

    if (size(inds) == 0) {{
        sim.setValue("meanFitness", NAN);
        return;
    }}

    age_survival = ageSurvivalForIndividuals(inds);
    inds.fitnessScaling = age_survival;

    denom = pop.individualCount * mean(age_survival);

    if (denom > 0.0)
        pop.fitnessScaling = targetN / denom;
    else
        pop.fitnessScaling = 0.0;

    sim.recalculateFitness();

    fit = pop.cachedFitness(NULL);
    if (size(fit) > 0)
        sim.setValue("meanFitness", mean(fit));
    else
        sim.setValue("meanFitness", NAN);
}}

function (void)enforcePopulationCeiling(object<Subpopulation>$ pop, integer$ maxN)
{{
    if (pop.individualCount > maxN) {{
        nRemove = pop.individualCount - maxN;
        sim.killIndividuals(pop.sampleIndividuals(nRemove));
    }}
}}

function (void)applyFutureCeilingFitnessScaling(object<Subpopulation>$ pop, numeric$ Kceiling)
{{
    inds = pop.individuals;

    if (size(inds) == 0) {{
        sim.setValue("meanFitness", NAN);
        return;
    }}

    age_survival = ageSurvivalForIndividuals(inds);
    inds.fitnessScaling = age_survival;

    // Future phase: upper carrying-capacity ceiling only.
    // If N is below Kceiling, there is no density-dependent rescue.
    // If N is above Kceiling, survival is down-scaled, but never up-scaled.
    if (pop.individualCount > Kceiling) {{
        denom = pop.individualCount * mean(age_survival);

        if (denom > 0.0)
            pop.fitnessScaling = min(c(1.0, Kceiling / denom));
        else
            pop.fitnessScaling = 0.0;
    }} else {{
        pop.fitnessScaling = 1.0;
    }}

    sim.recalculateFitness();

    fit = pop.cachedFitness(NULL);
    if (size(fit) > 0)
        sim.setValue("meanFitness", mean(fit));
    else
        sim.setValue("meanFitness", NAN);
}}

function (float$)meanHetAcrossChromosomesForIndividuals(object<Individual> inds)
{{
    if (size(inds) == 0)
        return NAN;

    totalHet = 0.0;
    totalLength = 0;

    for (chr in sim.chromosomes) {{
        haps = c();

        for (ind in inds)
            haps = c(haps, ind.haplosomesForChromosomes(chr, includeNulls=F));

        if (size(haps) > 1) {{
            totalHet = totalHet + calcHeterozygosity(haps) * chr.length;
            totalLength = totalLength + chr.length;
        }}
    }}

    if (totalLength == 0)
        return NAN;

    return totalHet / totalLength;
}}


function (float)mutationClassCountsForIndividuals(object<Individual> inds)
{{
    if (size(inds) == 0)
        return c(NAN, NAN, NAN, NAN, NAN, NAN, NAN);

    neutralCounts = c();
    weakDelCounts = c();
    moderateDelCounts = c();
    strongDelCounts = c();
    veryStrongDelCounts = c();
    lethalDelCounts = c();
    totalDelCounts = c();

    for (ind in inds) {{
        muts = c();

        for (chr in sim.chromosomes) {{
            haps = ind.haplosomesForChromosomes(chr, includeNulls=F);
            muts = c(muts, haps.mutations);
        }}

        if (size(muts) == 0) {{
            neutralCounts = c(neutralCounts, 0);
            weakDelCounts = c(weakDelCounts, 0);
            moderateDelCounts = c(moderateDelCounts, 0);
            strongDelCounts = c(strongDelCounts, 0);
            veryStrongDelCounts = c(veryStrongDelCounts, 0);
            lethalDelCounts = c(lethalDelCounts, 0);
            totalDelCounts = c(totalDelCounts, 0);
            next;
        }}

        nNeutral = sum(muts.mutationType == m1);
        nWeak = sum(muts.mutationType == m2);
        nModerate = sum(muts.mutationType == m3);
        nStrong = sum(muts.mutationType == m4);
        nVeryStrong = sum(muts.mutationType == m5);
        nLethal = sum(muts.mutationType == m6);
        nTotalDel = nWeak + nModerate + nStrong + nVeryStrong + nLethal;

        neutralCounts = c(neutralCounts, nNeutral);
        weakDelCounts = c(weakDelCounts, nWeak);
        moderateDelCounts = c(moderateDelCounts, nModerate);
        strongDelCounts = c(strongDelCounts, nStrong);
        veryStrongDelCounts = c(veryStrongDelCounts, nVeryStrong);
        lethalDelCounts = c(lethalDelCounts, nLethal);
        totalDelCounts = c(totalDelCounts, nTotalDel);
    }}

    return c(mean(neutralCounts),
             mean(weakDelCounts),
             mean(moderateDelCounts),
             mean(strongDelCounts),
             mean(veryStrongDelCounts),
             mean(lethalDelCounts),
             mean(totalDelCounts));
}}

function (void)writeHeader(void)
{{
    header = "gen\\tKceiling\\tpopSize\\tsampledN\\tmeanHet\\tFROH\\tmeanFitness\\tmeanNeutralMutCount\\tmeanWeakDelCount\\tmeanModerateDelCount\\tmeanStrongDelCount\\tmeanVeryStrongDelCount\\tmeanLethalDelCount\\tmeanTotalDelCount";
    writeFile(OUT_FILE, header, append=F);
}}

function (void)logStats(integer$ futureGen, object<Subpopulation>$ pop)
{{
    if (pop.individualCount == 0) {{
        line = paste(c(futureGen, FUTURE_K_CEILING, 0, 0, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN), sep="\\t");
        writeFile(OUT_FILE, line, append=T);
        return;
    }}

    // Recalculate fitness at the exact logging point.
    applyFutureCeilingFitnessScaling(pop, FUTURE_K_CEILING);

    nSample = asInteger(min(c(STATS_SAMPLE_SIZE, pop.individualCount)));
    sampledInds = pop.sampleIndividuals(nSample);

    meanHet = meanHetAcrossChromosomesForIndividuals(sampledInds);
    meanFroh = calcMeanFroh(sampledInds, minimumLength=ROH_CUTOFF_BP);

    sampleFitness = pop.cachedFitness(sampledInds.index);
    meanFitness = mean(sampleFitness);

    mutationCounts = mutationClassCountsForIndividuals(sampledInds);

    line = paste(c(futureGen,
                   FUTURE_K_CEILING,
                   pop.individualCount,
                   nSample,
                   meanHet,
                   meanFroh,
                   meanFitness,
                   mutationCounts), sep="\\t");

    writeFile(OUT_FILE, line, append=T);
}}

reproduction(NULL, "F")
{{
    if (individual.age < AGE_FIRST_REPRO)
        return;

    fit = subpop.cachedFitness(NULL)[individual.index];

    if (isNAN(fit))
        fit = 0.0;

    fit = max(c(0.0, fit));

    // Adult breeding probability is drawn for each adult female in each year
    // from a normal distribution with mean 0.8 and SD 0.1, truncated to [0, 1].
    baseReproProb = rnorm(1, ADULT_REPRO_PROB_MEAN, ADULT_REPRO_PROB_SD);
    baseReproProb = min(c(1.0, max(c(0.0, baseReproProb))));

    p_repro = baseReproProb * min(c(1.0, fit));

    if (runif(1) > p_repro)
        return;

    mate = subpop.sampleIndividuals(1, sex="M", minAge=AGE_FIRST_REPRO);

    if (size(mate) == 0)
        return;

    clutch = asInteger(round(rnorm(1, CLUTCH_MEAN, CLUTCH_SD)));
    clutch = max(c(0, clutch));

    if (clutch > 0) {{
        // Convert clutch size into recruited offspring entering the simulated population.
        // Each egg has probability CLUTCH_TO_ADULT_PROB (=0.39) of becoming a recruited individual.
        nRecruits = rbinom(1, clutch, CLUTCH_TO_ADULT_PROB);

        if (nRecruits > 0)
            subpop.addCrossed(individual, mate, count=nRecruits);
    }}
}}

survival()
{{
    if (runif(1) < sim.getValue("p_death"))
        return F;

    return NULL;
}}

1 early()
{{
    sim.addSubpop("p1", ANCESTRAL_N, sexRatio=0.5);
    p1.individuals.age = rdunif(p1.individualCount, min=0, max=MAX_LIFESPAN - 1);

    sim.setValue("p_death", 0.0);
    sim.setValue("meanFitness", NAN);
    sim.setValue("phase", "burnin");
    sim.setValue("futureGen", -1);
}}

2:{end_tick} early()
{{
    if ((sim.cycle <= BURNIN_GENERATIONS) & ((sim.cycle % BURNIN_PROGRESS_INTERVAL) == 0))
        catn("burn-in progress: tick=" + sim.cycle + ", N_p1=" + p1.individualCount);

    if (sim.getValue("phase") == "future")
        sim.setValue("p_death", rbeta(1, DISASTER_ALPHA, DISASTER_BETA));
    else
        sim.setValue("p_death", 0.0);

    if (exists("p1"))
        applyRegulatedFitnessScaling(p1, ANCESTRAL_N);

    if (exists("p3")) {{
        if (sim.getValue("phase") == "current")
            applyRegulatedFitnessScaling(p3, CURRENT_N);
        else if (sim.getValue("phase") == "future")
            applyFutureCeilingFitnessScaling(p3, FUTURE_K_CEILING);
    }}
	if ((sim.cycle % BURNIN_PROGRESS_INTERVAL) == 0) {{
        nP1 = -1;
        nP3 = -1;

        if (exists("p1"))
            nP1 = p1.individualCount;

        if (exists("p3"))
            nP3 = p3.individualCount;

        catn("tick=" + sim.cycle +
             ", phase=" + sim.getValue("phase") +
             ", N_p1=" + nP1 +
             ", N_p3=" + nP3 +
             ", mutations=" + size(sim.mutations) +
             ", substitutions=" + size(sim.substitutions));
    }}
}}

2:{burnin_generations - 1} late()
{{
    if (exists("p1") & (sim.getValue("phase") == "burnin"))
        enforcePopulationCeiling(p1, ANCESTRAL_N);
}}

{burnin_generations} late()
{{
    // Hard ceiling for p1 before drawing founders.
    enforcePopulationCeiling(p1, ANCESTRAL_N);

    nFounders = min(c(CURRENT_N, p1.individualCount));

    sim.addSubpop("p3", 0, sexRatio=0.5);

    if (nFounders > 0)
        p3.takeMigrants(p1.sampleIndividuals(nFounders));

    if (p3.individualCount > 0)
        p3.individuals.age = rdunif(p3.individualCount, min=0, max=MAX_LIFESPAN - 1);

    // No future supplementation is used in this version.
    p1.removeSubpopulation();

    sim.setValue("phase", "current");
    sim.setValue("currentGen", 0);

    catn("CURRENT START: " + SPECIES_NAME + " " + SCENARIO +
         ", p3 N=" + p3.individualCount +
         ", future K ceiling=" + FUTURE_K_CEILING);
}}

{burnin_generations + 1}:{burnin_generations + CURRENT_PHASE_GENERATIONS - 1} late()
{{
    if (exists("p3") & (sim.getValue("phase") == "current"))
        enforcePopulationCeiling(p3, CURRENT_N);
}}

{burnin_generations + CURRENT_PHASE_GENERATIONS} late()
{{
    // Hard ceiling for p3 at the end of the current phase.
    enforcePopulationCeiling(p3, CURRENT_N);

    // Future starts from the current p3 size.
    // FUTURE_K_CEILING is only an upper carrying-capacity ceiling,
    // not an initial population size and not a lower-bound target.
    sim.setValue("phase", "future");
    sim.setValue("futureGen", 0);

    writeHeader();

    applyFutureCeilingFitnessScaling(p3, FUTURE_K_CEILING);
    logStats(0, p3);

    catn("FUTURE START: " + SPECIES_NAME + " " + SCENARIO +
         ", actual initial N=" + p3.individualCount +
         ", future K ceiling=" + FUTURE_K_CEILING +
         ", no lower-bound rescue");
}}

{future_start_tick}:{end_tick} late()
{{
    if (sim.getValue("phase") != "future")
        return;

    // Hard ceiling in the future phase:
    // if p3 exceeds FUTURE_K_CEILING after reproduction/survival,
    // randomly remove excess individuals.
    enforcePopulationCeiling(p3, FUTURE_K_CEILING);

    futureGen = sim.getValue("futureGen") + 1;
    sim.setValue("futureGen", futureGen);

    if ((futureGen % LOG_INTERVAL) == 0)
        logStats(futureGen, p3);

    if (p3.individualCount == 0) {{
        if ((futureGen % LOG_INTERVAL) != 0)
            logStats(futureGen, p3);

        catn("EXTINCT: " + SPECIES_NAME + " " + SCENARIO +
             " at future generation " + futureGen);

        sim.simulationFinished();
    }} else if (futureGen >= FUTURE_MAX_GENERATIONS) {{
        catn("FINISHED: " + SPECIES_NAME + " " + SCENARIO +
             " reached " + FUTURE_MAX_GENERATIONS + " future generations");

        sim.simulationFinished();
    }}
}}
'''
    return textwrap.dedent(template).strip() + "\n"


def main() -> None:
    MODEL_DIR.mkdir(exist_ok=True)
    RESULT_DIR.mkdir(exist_ok=True)

    model_paths: List[Path] = []

    for sp in SPECIES_TABLE:
        safe_species = sanitize_name(sp.species)

        for scenario_name, multiplier in SCENARIOS.items():
            model_path = MODEL_DIR / f"{safe_species}.{scenario_name}.slim"
            model_path.write_text(
                build_slim_script(sp, scenario_name, multiplier),
                encoding="utf-8",
            )
            model_paths.append(model_path)

    run_lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        f"mkdir -p {RESULT_DIR.as_posix()}",
        "",
    ]

    for model_path in model_paths:
        run_lines.append(f"slim {model_path.as_posix()}")

    RUN_SCRIPT.write_text("\n".join(run_lines) + "\n", encoding="utf-8")
    RUN_SCRIPT.chmod(0o755)

    print(f"Generated {len(model_paths)} SLiM model(s) in {MODEL_DIR}/")
    print(f"Generated runner: {RUN_SCRIPT}")
    print(f"Results will be written to {RESULT_DIR}/")
    print("Burn-in phase is regulated around ancestral_N.")
    print("Current phase is regulated around current_N.")
    print("Future scenarios set only an upper carrying-capacity ceiling.")
    print("Future phase starts from the current p3 size; no one-off resizing is used.")
    print("No lower-bound rescue is applied when future N is below the ceiling.")
    print(f"Current phase generations: {CURRENT_PHASE_GENERATIONS}")
    print(f"Future output interval: every {LOG_INTERVAL} generations")
    print("Adult breeding probability: N(0.8, 0.1), truncated to [0, 1].")
    print("Each egg has probability 0.39 of becoming a recruited offspring.")
    print("Neutral mutations are simulated; neutral:deleterious = 1:2.31.")
    print("Age-specific annual survival is applied through LIFE_TABLE and individual fitnessScaling.")


if __name__ == "__main__":
    main()
