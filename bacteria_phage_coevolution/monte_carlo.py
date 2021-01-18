from collections import namedtuple
from math import log
from random import random


GenState = namedtuple(
    "GenState",
    [
        "uninfected_bacteria",
        "phages",
        "infected_bacteria",
        "carry_over",
        "gen_count",
        "past_total_uninfected",
    ],
)
MCParameters = namedtuple(
    "MCParameters",
    ["adsorption_rate", "burst_size", "latent_period"],
)


def calculate_next_gen(CurrentState, SimulationParameters):
    (
        uninfected_bacteria,
        phages,
        infected_bacteria,
        carry_over,
        gen_count,
        past_total_uninfected,
    ) = CurrentState

    second_last_uninfected, last_uninfected = past_total_uninfected
    total_uninfected = last_uninfected

    dN = last_uninfected - second_last_uninfected
    duplication_rate = log(dN) / 10

    for _ in range(uninfected_bacteria):
        random_1 = random()
        random_2 = random()

        if random_1 >= SimulationParameters.adsorption_rate:
            if random_2 >= duplication_rate:
                continue
            elif random_2 < duplication_rate:
                uninfected_bacteria += 1
                total_uninfected += 1
        elif random_1 < SimulationParameters.adsorption_rate:
            uninfected_bacteria -= 1
            infected_bacteria += 1
            carry_over[gen_count + SimulationParameters.latent_period][
                "phage_count"
            ] += SimulationParameters.burst_size
            carry_over[gen_count + SimulationParameters.latent_period][
                "infected_bacteria_count"
            ] -= 1

    past_total_uninfected = (last_uninfected, total_uninfected)

    NextGen = GenState(
        uninfected_bacteria,
        phages,
        infected_bacteria,
        carry_over,
        gen_count,
        past_total_uninfected,
    )
    return NextGen


def apply_carry_over(NextGen):
    phages = NextGen.phages
    infected_bacteria = NextGen.infected_bacteria
    gen_count = NextGen.gen_count

    phages += NextGen.carry_over[gen_count]["phage_count"]
    infected_bacteria += NextGen.carry_over[gen_count]["infected_bacteria_count"]

    gen_count += 1
    UpdatedNextGen = GenState(
        NextGen.uninfected_bacteria,
        phages,
        infected_bacteria,
        NextGen.carry_over,
        gen_count,
        NextGen.past_total_uninfected,
    )
    return UpdatedNextGen


def monte_carlo(n):
    SimulationParameters = MCParameters(0.5, 10, 2)

    carry_over = [
        {"phage_count": 0, "infected_bacteria_count": 0}
        for _ in range(n + SimulationParameters.latent_period)
    ]
    GS = GenState(1000, 500, 0, carry_over, 0, (1000, 1020))

    for _ in range(n):
        GS = apply_carry_over(calculate_next_gen(GS, SimulationParameters))
    return GS


print(monte_carlo(10))
