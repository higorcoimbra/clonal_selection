from random import random, uniform
from math import sin, sqrt, exp
#import matplotlib.pyplot as plt
from operator import attrgetter

offset = 1500

# function domain
min_domain = -512
max_domain = 512

class Antibody:
	def __init__(self, x1, x2, affinity):
		self.x1 = x1
		self.x2 = x2
		self.affinity = affinity


# define population list
def generate_population(pop_size):
	population = []
	for i in range(0,pop_size):
		antibody = Antibody(uniform(min_domain,max_domain),uniform(min_domain,max_domain),0)
		population.append(antibody)
	return population

# generate affinity value for every antibody in population
def affinity(population):
	global offset
	for c in population:
		c.affinity = -(c.x2+47)*sin(sqrt(abs(c.x2+c.x1/2+47)))-c.x1*sin(sqrt(abs(c.x1-(c.x2+47))))
		c.affinity = -c.affinity+offset
	return population	

def cloning(population,num_clones):
	# generating num_clones clones for every antibody available
	clones = []
	for antibody in population:
		for i in range(0,num_clones):
			clones.append(antibody)
	return clones

def mutation(population,ro):
	global max_domain
	global min_domain
	max_affinity = max(population,key=attrgetter('affinity')).affinity
	for c in population:
		d_star = c.affinity/max_affinity
		mutation_rate = exp(-ro*d_star)
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

def selection(clone_population, num_clones):
	
	cluster = []
	population = []
	i = 0
	while(i < len(clone_population)):
		cluster = clone_population[i:i+num_clones]
		selected_clone = max(cluster,key=attrgetter('affinity'))
		population.append(selected_clone)
		cluster.clear()
		i = i+num_clones

	return population

def print_population(population):
	for c in population:
		print(str(c.x1)+" "+str(c.x2)+" affinity: "+str(c.affinity))

pop_size = 500
ro = 0.2
num_clones = 5
population = generate_population(pop_size)
iterations = 50

# multi modal optimization process
for i in range(0,iterations):
	population = affinity(population)
	# general process of clonal algorithm
	clone_population = cloning(population,num_clones)
	clone_population = mutation(clone_population,ro)
	population = selection(clone_population,num_clones)

print_population(population)

# plot for statistical analysis

