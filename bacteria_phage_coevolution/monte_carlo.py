from collections import namedtuple
from math import e, exp
from random import random


GenState = namedtuple(
    "GenState",
    [
        "uninfected_bacteria",
        "phages",
        "infected_bacteria",
        "carry_over",
        "gen_count",
    ],
)
MCParameters = namedtuple(
    "MCParameters",
    ["growth_rate", "carrying_capacity", "burst_size", "latent_period"],
)


def calculate_next_gen(CurrentState, SimulationParameters):
    (
        uninfected_bacteria,
        phages,
        infected_bacteria,
        carry_over,
        _,
    ) = CurrentState

    increased_uninfected_bacteria = (
        SimulationParameters.growth_rate
        * uninfected_bacteria
        * (1 - uninfected_bacteria / SimulationParameters.carrying_capacity)
    )
    dN = increased_uninfected_bacteria - uninfected_bacteria
    duplication_rate = dN / uninfected_bacteria

    moi = phages / uninfected_bacteria
    adsorption_rate = moi * exp(-moi)

    carry_over_index = CurrentState.gen_count + SimulationParameters.latent_period
    carry_over_size = len(carry_over)

    carry_over[carry_over_index % carry_over_size] = {
        "phage_count": 0,
        "infected_bacteria_count": 0,
        "uninfected_bacteria_count": 0,
    }

    for _ in range(uninfected_bacteria):
        random_1 = random()
        random_2 = random()
        if random_1 >= adsorption_rate:
            if random_2 >= duplication_rate:
                continue
            elif random_2 < duplication_rate:
                uninfected_bacteria += 1
        elif random_1 < adsorption_rate:
            uninfected_bacteria -= 1
            infected_bacteria += 1
            phages -= 1

            carry_over[carry_over_index % carry_over_size][
                "phage_count"
            ] += SimulationParameters.burst_size
            carry_over[carry_over_index % carry_over_size][
                "infected_bacteria_count"
            ] -= 1

            a = e/2
            secondary_killing = int(
                uninfected_bacteria
                * exp(
                    -CurrentState.gen_count / SimulationParameters.latent_period * a
                )
            )
            carry_over[carry_over_index % carry_over_size][
                "uninfected_bacteria_count"
            ] -= secondary_killing

    carry_over[carry_over_index % carry_over_size]["phage_count"] //= 20

    NextGen = GenState(
        uninfected_bacteria,
        phages,
        infected_bacteria,
        carry_over,
        CurrentState.gen_count,
    )
    return NextGen


def apply_carry_over(NextGen):
    phages = NextGen.phages
    infected_bacteria = NextGen.infected_bacteria
    uninfected_bacteria = NextGen.uninfected_bacteria
    gen_count = NextGen.gen_count

    carry_over_size = len(NextGen.carry_over)

    phages += NextGen.carry_over[gen_count % carry_over_size]["phage_count"]
    infected_bacteria += NextGen.carry_over[gen_count % carry_over_size][
        "infected_bacteria_count"
    ]

    secondary_killing_rate = abs(NextGen.carry_over[gen_count % carry_over_size]["uninfected_bacteria_count"] / uninfected_bacteria)
    print(secondary_killing_rate)
    for _ in range(uninfected_bacteria):
        r = random()
        if r < secondary_killing_rate:
            uninfected_bacteria -= 1

    gen_count += 1
    UpdatedNextGen = GenState(
        NextGen.uninfected_bacteria,
        phages,
        infected_bacteria,
        NextGen.carry_over,
        gen_count,
    )
    return UpdatedNextGen


def monte_carlo(n):
    SimulationParameters = MCParameters(3, 5000, 150, 3)

    carry_over = [
        {"phage_count": 0, "infected_bacteria_count": 0, "uninfected_bacteria_count": 0}
        for _ in range(SimulationParameters.latent_period + 1)
    ]
    GS = GenState(1000, 1000, 0, carry_over, 0)

    for _ in range(n):
        GS = apply_carry_over(calculate_next_gen(GS, SimulationParameters))
    return GS


FinalState = monte_carlo(100)
print(f"Uninfected bacteria count: {FinalState.uninfected_bacteria}")
print(f"Phages count: {FinalState.phages}")
print(f"Infected bacteria: {FinalState.infected_bacteria}")
