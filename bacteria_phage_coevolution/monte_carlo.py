from collections import namedtuple
from random import random


GenState = namedtuple(
    "GenState", ["uninfected_bacteria", "phages", "infected_bacteria"]
)
MCParameters = namedtuple(
    "MCParameters", ["adsorption_rate", "duplication_rate", "burst_size"]
)


def calculate_next_gen(CurrentState, SimulationParameters):
    uninfected_bacteria, phages, infected_bacteria = CurrentState
    for _ in range(uninfected_bacteria):
        random_1 = random()
        random_2 = random()

        if random_1 >= SimulationParameters.adsorption_rate:
            if random_2 >= SimulationParameters.duplication_rate:
                continue
            elif random_2 < SimulationParameters.duplication_rate:
                uninfected_bacteria += 1
        elif random_1 < SimulationParameters.adsorption_rate:
            uninfected_bacteria -= 1
            infected_bacteria += 1
            phages += SimulationParameters.burst_size

    NextGen = GenState(uninfected_bacteria, phages, infected_bacteria)
    return NextGen


def monte_carlo(n):
    SimulationParameters = MCParameters(0.5, 0.3, 10)
    GS = GenState(1000, 500, 0)
    for _ in range(n):
        GS = calculate_next_gen(GS, SimulationParameters)
    return GS
