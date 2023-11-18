import random
from deap import base
from deap import creator
from deap import tools
from deap import algorithms
import pandas as pd
from rdkit import rdBase, Chem, DataStructs
from numpy import dot 
from numpy.linalg import norm 
import itertools
import pubchempy as pcp
from rdkit import rdBase, Chem, DataStructs
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.AtomPairs import Pairs, Torsions
import numpy
from rdkit import RDLogger 
import time

trys = 25

now = time.ctime()
cnvtime = time.strptime(now)
print(time.strftime("%Y/%m/%d %H:%M", cnvtime), flush=True)

RDLogger.DisableLog('rdApp.*') 

df = pd.read_csv("231103_use.tsv", sep="\t")

creator.create("Fitness", base.Fitness, weights=(-1,-2,6))
creator.create("Individual", set, fitness=creator.Fitness)

IND_INIT_SIZE = 10
toolbox = base.Toolbox()
toolbox.register("attr_item", random.randint, 0, len(df)-1)
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_item, IND_INIT_SIZE)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

def cxSet(ind1, ind2):
    ind_1 = list(ind1)
    ind_2 = list(ind2)
    a = random.randint(0,IND_INIT_SIZE-1)
    ind1 = ind_1[:a] + ind_2[a:]
    ind2 = ind_1[a:] + ind_2[:a]
    ind1 = creator.Individual(ind1)
    ind2 = creator.Individual(ind2)  
    return ind1, ind2

def not_dup(a, b, list):
    k = 0
    while k < 1:
        s = random.randint(a,b)
        if s not in list:
            k = 100
    return s

def mutSet(individual):
    for i in range(len(individual)):
        ind_list = list(individual)
        if len(ind_list) != IND_INIT_SIZE:
            return 1
        if random.random() < 1:
            n = random.randint(0,IND_INIT_SIZE-1)
            ind_list[n] = not_dup(0, len(df)-1, ind_list)
            individual = creator.Individual(ind_list)
        
    return individual

def multi_oo(individual):
    if len(individual) != IND_INIT_SIZE:
        return 100000, 100000, -100000

    total_tanimoto = 0
    total_cosine = 0
    total_toxicity = 0
    total_MW = 0

    t = list(itertools.combinations(individual,2))

    for i in range(len(t)):
        if t[i][0] == t[i][1]:
            print("!", flush=True)
            return 100000, 100000, -100000

    for i in range(len(t)):
        mol1 = Chem.MolFromSmiles(df.iloc[t[i][0],5])
        fp_1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius=2, nBits=2048)
        mol2 = Chem.MolFromSmiles(df.iloc[t[i][1],5])
        fp_2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius=2, nBits=2048)
        tanimoto = DataStructs.TanimotoSimilarity(fp_1, fp_2)

        a = numpy.array([df.iloc[t[i][0],7],df.iloc[t[i][0],8]], dtype=float)
        b = numpy.array([df.iloc[t[i][1],7],df.iloc[t[i][1],8]], dtype=float)
        cos = dot(a, b) / (norm(a) * norm(b))
        
        c = df.iloc[t[i][0], 4]
        d = df.iloc[t[i][1], 4]
        MW = abs(c - d)

        tox = (df.iloc[t[i][0],1] - df.iloc[t[i][1],1]) ** 2

        total_tanimoto += tanimoto
        total_cosine += cos
        total_toxicity += tox
        #total_MW += MW
    
    ind_list = list(individual)
    sev_0 = df.iloc[ind_list[0], 1]
    sev_1 = df.iloc[ind_list[1], 1]
    sev_2 = df.iloc[ind_list[2], 1]
    sev_3 = df.iloc[ind_list[3], 1]
    sev_4 = df.iloc[ind_list[4], 1]
    sev_5 = df.iloc[ind_list[5], 1]
    sev_6 = df.iloc[ind_list[6], 1]
    sev_7 = df.iloc[ind_list[7], 1]
    sev_8 = df.iloc[ind_list[8], 1]
    sev_9 = df.iloc[ind_list[9], 1]
    sev_set = {sev_0, sev_1, sev_2, sev_3, sev_4, sev_5, sev_6, sev_7, sev_8, sev_9}
    print(sev_set)
    total_toxicity += (len(sev_set) * 100) 
    #if len(sev_set) == 9:
    #    total_toxicity += 10000
 
    return total_tanimoto, total_cosine, total_toxicity #total_MW

toolbox.register("evaluate", multi_oo)
toolbox.register("mate", cxSet)
#toolbox.register("mate", cxTwoPointCopy)
toolbox.register("mutate", mutSet)
toolbox.register("select", tools.selNSGA2)

def main():
    random.seed(42)
    N_GEN = 100
    POP_SIZE = 100
    CX_PB = 0.5
    MUT_PB = 0.2
    tsv = []
    tsv_2 = []
    tsv_3 = []

    pop = toolbox.population(n=POP_SIZE)

    fitnesses = list(map(toolbox.evaluate, pop))
    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit
    #print("  Evaluated %i individuals" % len(pop))

    fits = [ind.fitness.values[0] for ind in pop]

    g = 0
    while g < N_GEN:
        g = g + 1
        if g % 30 == 0:
            print("!", flush=True)
        #print("-- Generation %i --" % g)
        offspring = toolbox.select(pop, 100)
        offspring = list(map(toolbox.clone, offspring))
        #print(len(offspring))

        for child1, child2 in zip(offspring[::2], offspring[1::2]):
            if random.random() < CX_PB:
                child1, child2 = toolbox.mate(child1, child2)
                del child1.fitness.values
                del child2.fitness.values
            
            if len(child1) == IND_INIT_SIZE: 
                offspring.append(child1)
            if len(child2) == IND_INIT_SIZE:
                offspring.append(child2)
        #print(offspring)

        for mutant in offspring:
            if random.random() < MUT_PB:
                mutant = toolbox.mutate(mutant)
                if mutant == 1:
                    continue
                offspring.append(mutant)

        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        #print("  Evaluated %i individuals" % len(invalid_ind))
        
        #print(f"BEST: {tools.selBest(offspring,1)}")
        tsv_2.append(tools.selBest(offspring,5))

        pop[:] = offspring
        #print(pop)

        col = []
        # 新世代の全個体の適応度の抽出
        fits_0 = [ind.fitness.values[0] for ind in pop]

        # 適応度の統計情報の出力
        length = len(pop)
        mean = sum(fits_0) / length
        sum2 = sum(x*x for x in fits_0)
        std = abs(sum2 / length - mean**2)**0.5
        col.append(mean)
        col.append(std)
        #print("===TANIMOTO===")
        #print("  Min %s" % min(fits_0))
        #print("  Max %s" % max(fits_0))
        #print("  Avg %s" % mean)
        #print("  Std %s" % std)

        # 新世代の全個体の適応度の抽出
        fits_1 = [ind.fitness.values[1] for ind in pop]

        # 適応度の統計情報の出力
        length = len(pop)
        mean = sum(fits_1) / length
        sum2 = sum(x*x for x in fits_1)
        std = abs(sum2 / length - mean**2)**0.5
        col.append(mean)
        col.append(std)
        #print("===BUSSEI===")
        #print("  Min %s" % min(fits_1))
        #print("  Max %s" % max(fits_1))
        #print("  Avg %s" % mean)
        #print("  Std %s" % std)

        # 新世代の全個体の適応度の抽出
        fits_2 = [ind.fitness.values[2] for ind in pop]

        # 適応度の統計情報の出力
        length = len(pop)
        mean = sum(fits_2) / length
        sum2 = sum(x*x for x in fits_2)
        std = abs(sum2 / length - mean**2)**0.5
        col.append(mean)
        col.append(std)
        #print("===TOXICITY===")
        #print("  Min %s" % min(fits_2))
        #print("  Max %s" % max(fits_2))
        #print("  Avg %s" % mean)
        #print("  Std %s" % std)

        # 新世代の全個体の適応度の抽出
        #fits_3 = [ind.fitness.values[3] for ind in pop]

        # 適応度の統計情報の出力
        #length = len(pop)
        #mean = sum(fits_3) / length
        ##sum2 = sum(x*x for x in fits_3)
        #std = abs(sum2 / length - mean**2)**0.5
        #col.append(mean)
        #col.append(std)

        tsv.append(col)
      
        tsv_3.append(pop)

    

    # 最良個体の抽出
    best_ind = tools.selNSGA2(pop, 1)[0]
    print("Best individual is %s, %s" % (best_ind, best_ind.fitness.values))
    
    pd.DataFrame(tsv).to_csv(f"231102_result_{trys}_score.tsv", sep="\t")
    pd.DataFrame(tsv_2).to_csv(f"231102_result_{trys}_top_5.tsv", sep="\t")
    pd.DataFrame(tsv_3).to_csv(f"231102_result_{trys}_pop.tsv", sep="\t")

if __name__ == "__main__":
    main()

now = time.ctime()
cnvtime = time.strptime(now)
print(time.strftime("%Y/%m/%d %H:%M", cnvtime), flush=True)