from random import random, uniform
from math import sin, sqrt
import matplotlib.pyplot as plt

# function domain
min_domain = -512
max_domain = 512

class Cromossome:
	def __init__(self, x1, x2,fitness):
		self.x1 = x1
		self.x2 = x2
		self.fitness = fitness


# define population list
def generate_population(pop_size):
	population = []
	for i in range(0,pop_size):
		cromossome = Cromossome(uniform(min_domain,max_domain),uniform(min_domain,max_domain),0)
		population.append(cromossome)
	return population

# generate fitness value for every cromossome in population
def fitness(population):
	max_fitness = 0
	min_fitness = 2500
	mid_fitness = 0
	maior = Cromossome(0, 0, 0)
	for c in population:
		c.fitness = -(c.x2+47)*sin(sqrt(abs(c.x2+c.x1/2+47)))-c.x1*sin(sqrt(abs(c.x1-(c.x2+47))))
		c.fitness = -c.fitness+offset
		if (c.fitness > max_fitness):
			max_fitness = c.fitness
			maior = Cromossome(c.x1, c.x2, c.fitness)
		min_fitness = min(c.fitness, min_fitness)
		mid_fitness += c.fitness 

	mid_fitness /= len(population)
	return [population, max_fitness, min_fitness, mid_fitness, maior]

# selection process - roulette wheel
def selection(population,pop_size):
	roullete = []
	value = 0
	# roullete ranges
	roullete.append(0)
	for c in population:
		value += c.fitness
		roullete.append(value)

	# selection with random number
	selected_parents = []
	for i in range(0,pop_size):
		rand = uniform(0.0,max(roullete))
		for i in range(1,len(roullete)):
			if(rand > roullete[i-1] and rand <= roullete[i]):
				selected_parents.append(i-1)
	return selected_parents

def normalize(c):
	global min_domain
	global max_domain
	if(c.x1 < min_domain):
		c.x1 = min_domain
	if(c.x1 > max_domain):
		c.x1 = max_domain
	if(c.x2 < min_domain):
		c.x2 = min_domain
	if(c.x2 > max_domain):
		c.x2 = max_domain
	return c
	

def crossover(population, selected_parents, crossover_rate,elite,elite_bool):
	sons = []
	i = 0
	while i < len(selected_parents):
		alpha = random()
		parents = [population[selected_parents[i]], population[selected_parents[i+1]]]
		if(i == 0 and elite_bool):
			son1 = elite
			son2 = elite
		elif(random() < crossover_rate):
			son1 = Cromossome(0,0,0)
			son2 = Cromossome(0,0,0)
			son1.x1 = alpha*parents[0].x1+(1-alpha)*parents[1].x1
			son1.x2 = alpha*parents[0].x2+(1-alpha)*parents[1].x2
			son2.x1 = alpha*parents[1].x1+(1-alpha)*parents[0].x1
			son2.x2 = alpha*parents[1].x2+(1-alpha)*parents[0].x2
		else:
			son1 = population[selected_parents[i]]
			son2 = population[selected_parents[i+1]]
		son1 = normalize(son1)
		son2 = normalize(son2)
		sons.append(son1)
		sons.append(son2)
		i = i+2

	return sons

def mutacao(population, mutation_rate):
	global max_domain
	global min_domain
	for c in population:
		if(random() < mutation_rate):
			mutated_value = c.x1+uniform(-2.0,2.0)
			if(not(mutated_value < min_domain)):
				c.x1 = mutated_value
			else:
				c.x1 = min_domain
			if(not(mutated_value > max_domain)):
				c.x1 = mutated_value
			else:
				c.x1 = max_domain

		if(random() < mutation_rate):
			mutated_value = c.x2+uniform(-2.0,2.0)
			if(not(mutated_value < min_domain)):
				c.x2 = mutated_value
			else:
				c.x2 = min_domain
			if(not(mutated_value > max_domain)):
				c.x2 = mutated_value
			else:
				c.x2 = max_domain
	return population

def print_population(population):
	for c in population:
		print(str(c.x1)+" "+str(c.x2))

pop_size = 500
crossover_rate = 0.5
mutation_rate = 0.13
population = generate_population(pop_size)
generations = 200
max_fitness_list = []
mid_fitness_list = []
min_fitness_list = []

# evolution process
for i in range(0,generations):
	fitness_atribs = fitness(population)
	population = fitness_atribs[0]

	# array operations for statistical analysis
	max_fitness_list.append(fitness_atribs[1])
	min_fitness_list.append(fitness_atribs[2])
	mid_fitness_list.append(fitness_atribs[3])
	
	# general process of ga
	selected_parents = selection(population,pop_size)
	population = crossover(population, selected_parents, crossover_rate,fitness_atribs[4],True)
	population = mutacao(population,mutation_rate)

# individuo com maior fitness

print(fitness_atribs[4].x1)
print(fitness_atribs[4].x2)
print(fitness_atribs[4].fitness)

# plot for statistical analysis
plt.plot(max_fitness_list, label = 'Máximo fitness')
plt.plot(mid_fitness_list, label = 'Fitness médio')
plt.plot(min_fitness_list, label = 'Mínimo fitness')
plt.title("Evolução sem elitismo")
plt.xlabel("Gerações (n)")
plt.ylabel("Valor da função objetivo com operação -f(x)+1500")
plt.legend(loc = 2)
plt.show()

