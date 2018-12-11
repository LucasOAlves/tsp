import numpy as np
import copy
from sys import argv
from random import sample
from random import shuffle
from scipy.spatial import distance
import matplotlib.pyplot as plt

'''
Constants
'''

# Taxa de crossover
CONST_CROSS_RATE = 98

# Taxa de mutação
CONST_MUTATION_RATE = 4

# Tamanho da população
CONST_POPULATION_SIZE = 100

# Quantidade de gerações
CONS_GENERATIONS = 1000

# Quantidade de cidades
global CONST_CITIES
CONST_CITIES = 10

# Caso gerar aleatório limite dos eixos
CONST_LIMIT = 20

# Vetor com as cidades
global CITIES
CITIES = [sample(range((CONST_LIMIT)*-1,(CONST_LIMIT)),2) for x in range(CONST_CITIES)]

# Classe cromossomo
class Cromossomo:
	def __init__(self, cr=[]):
		# Verifica se a entrada e um cromossomo valido ou nao. 
		# Se for um vazio cria um novo cromossomo
		self.cro = self.make_initial_cro() if cr == [] else cr

		# Calcula o fitnes para o cromosso
		self.fit = self.calc_fitnes()

	# Funcao para criar um cromossomo aleatorio
	def make_initial_cro(self):

		# Cria um cromosso aleatorio com a quantidade de genes especificada
		city = copy.copy(CITIES)

		# Retorna este cromossomo criado
		shuffle(city)
		return city

	# Funcao para o calculo do fitnes
	def calc_fitnes(self):
		dist = 0 
		#concatena todos os genes em um unico para fazer a conta do fitnes
		for i in range(CONST_CITIES-1):
			# Soma as distancias euclidinas
			dist += distance.euclidean(self.cro[i],self.cro[i+1])
		
		# Soma a volta ao primeiro ponto
		dist += distance.euclidean(self.cro[CONST_CITIES-1],self.cro[0])

		# Retorna distancia
		return dist # calcula o fitnes

# Mostra a população
def show_population(pop):
	pop = sorted(pop,key = lambda x:x.fit)
	print('\n')
	for i in range(len(pop)):
		print('Ind: '+str(i)+'\t\tCro: '+str(pop[i].cro)+'\t\t\tFit: '+str(pop[i].fit))

# Função pra realizar o torneio
def make_tournament(pop):

	# Inicia população intermediaria vazio
	pop_middle = []

	# Salva o melhor individuo da população
	best = get_best(pop)

	# Para cada individuo da população
	for i in range(len(pop)):
		
		# sorteia dois indices
		c1,c2=np.random.choice((len(pop)),2,replace=False) #select 2 ind to make tournament

		# verifica qual o melhor para esses dois indices
		cr_max = min([pop[c1],pop[c2]],key = lambda x:x.fit) #verify what is the best

		# adiciona na população intermediaria
		pop_middle.append(Cromossomo(cr = cr_max.cro))
	
	# adiciona o melhor da populaçao anterior no pior da poulação intermediaria
	pop_middle[pop_middle.index(get_worst(pop_middle))] = copy.deepcopy(best)

	# retorna a população
	return pop_middle

def make_ox(pop):
	# inicia a população vazio
	pop_middle = []

	# salva o melhor da população
	best = get_best(pop)

	# para cada dois individuos da população
	for i in range(int(len(pop)/2)):

		#seleciona dois cromossomos da populacao para fazer o cross-over
		c1,c2 = np.random.choice((len(pop)),2, replace=False)
		
		#sorteia a problabilidade de fazer o cross-over
		prob = np.random.randint(0,101) 

		#verifica se a probabilidade sortiada caiu para fazer o cross-over
		if(prob<=CONST_CROSS_RATE): 

			# seleciona limite de substring
			cut = sample(range(1, CONST_CITIES-1),2)
			cut.sort();
			pf1 = [None]*CONST_CITIES
			pf2 = [None]*CONST_CITIES

			pf1[cut[0]:cut[1]] = pop[c1].cro[cut[0]:cut[1]]
			pf2[cut[0]:cut[1]] = pop[c2].cro[cut[0]:cut[1]]
			# print('-------------------------------------------------')
			# print(pf1,'||',pf2)
			f1Aux = copy.deepcopy(pop[c2])
			f2Aux = copy.deepcopy(pop[c1])
			# print(f1Aux.cro,'||',f2Aux.cro) 

			for j in range(len(pf1)):
				if(pf1[j] in f2Aux.cro):
					index = f2Aux.cro.index(pf1[j])
					del f2Aux.cro[index]

				if(pf2[j] in f1Aux.cro):
					index2 = f1Aux.cro.index(pf2[j])
					del f1Aux.cro[index2]
			count1 = 0
			count2 = 0
			for j in range(len(pf1)):
				if(pf1[j] == None):
					pf1[j] = f2Aux.cro[count1]
					count1 = count1 + 1
				if(pf2[j] == None):
					pf2[j] = f1Aux.cro[count2]
					count2 = count2 + 1

			# print('-------------------------------------------------')
			# print(pf1,'||',pf2)
			# print(f1Aux.cro,'||',f2Aux.cro) 
			# print('-------------------------------------------------')
			# input()
			pop_middle.append(Cromossomo (cr = pf1))
			pop_middle.append(Cromossomo (cr = pf2))
		else:
			# Se nao caiu na probabilidade adiciona os pais na população intermediária
			pop_middle.append(Cromossomo(cr = pop[c1].cro))
			pop_middle.append(Cromossomo(cr = pop[c2].cro))

	# substitui o pior da população intermediaria pelo melhor da poulação normal
	pop_middle[pop_middle.index(get_worst(pop_middle))] = copy.deepcopy(best)

	return pop_middle


# Realiza o OBX
def make_obx(pop):

	# inicia a população vazio
	pop_middle = []

	# salva o melhor da população
	best = get_best(pop)

	# para cada dois individuos da população
	for i in range(int(len(pop)/2)):

		#seleciona dois cromossomos da populacao para fazer o cross-over
		c1,c2 = np.random.choice((len(pop)),2, replace=False)
		
		#sorteia a problabilidade de fazer o cross-over
		prob = np.random.randint(0,101) 

		
		#verifica se a probabilidade sortiada caiu para fazer o cross-over
		if(prob<=CONST_CROSS_RATE): 

			#Pontos de ordem para o obx
			
			cut = sample(range(1, CONST_CITIES-1),int(CONST_CITIES/2))
			f1 = []
			f2 = []

			# encontra as ordem de pai1 em pai2 e vice versa
			ord1 = find_ord(cut,pop[c1].cro,pop[c2].cro)
			ord2 = find_ord(cut,pop[c2].cro,pop[c1].cro)
			
			# inicia contador
			cont = 0

			# para cada cidade
			for j in range(CONST_CITIES):

				# verifica se a cidade foi sorteada
				if(j in cut):	

					# adiciona na ordem correta em filho 1
					f1.append(pop[c2].cro[ord1[cont]])

					# adiciona na ordem correta em filho 2
					f2.append(pop[c1].cro[ord2[cont]])

					# soma 1 no contador de cidades
					cont = cont + 1

				else:
					# não esta nas cidades sorteadas somente adicona na ordem em que aparecem nos pais
					f1.append(pop[c1].cro[j])
					f2.append(pop[c2].cro[j])

			# adiciona os dois novos filhos gerados na populaçao intermediaria
			pop_middle.append(Cromossomo (cr = f1))
			pop_middle.append(Cromossomo (cr = f2))

		# Else de nao cair na probabilidade
		else:
			
			# Se nao caiu na probabilidade adiciona os pais na população intermediária
			pop_middle.append(Cromossomo(cr = pop[c1].cro))
			pop_middle.append(Cromossomo(cr = pop[c2].cro))
	
	# substitui o pior da população intermediaria pelo melhor da poulação normal
	pop_middle[pop_middle.index(get_worst(pop_middle))] = copy.deepcopy(best)

	# retorna população intermediaria
	return pop_middle

# Entra os as ordens
def find_ord(cut, p1, p2):

	# inicia auxiliar vazio
	aux = []

	# para ponto de corte
	for i in range(len(cut)):
		# adiciona na ordem que aparece no outro pai
		aux.append(p2.index(p1[cut[i]]))	
	# ordena
	aux.sort()
	# retorna aux
	return aux

# Realiza mutação por permutação de pontos
def make_ord_element_mutation(pop):	
	# salva o melhor individuo
	best = min(pop,key = lambda x:x.fit)

	# para individuo da poulação
	for i in range(len(pop)):

		# Sorteia probabilidade de fazer a mutacao
		prob = np.random.randint(0,101) 

		# testa a probabildiade
		if(prob < CONST_MUTATION_RATE):

			# sorteia os dois pontos de troca
			c1,c2 = np.random.choice(range(0,CONST_CITIES),2, replace=False)

			# salva o primeiro ponto em uma variavel aleatória
			aux = pop[i].cro[c1]

			# coloca no primeiro ponto o valor do segundo ponto
			pop[i].cro[c1] = pop[i].cro[c2]

			# coloca no segundo ponto ovalor salvo
			pop[i].cro[c2] = aux
	
	# troca pelo melhor da população inicial
	pop[pop.index(get_worst(pop))] = copy.deepcopy(best)

	# retorna pop
	return pop

# Realiza a mutação por inversão de sublista
def make_inversion_sublist_mutation(pop):

	# salva o melhor
	best = get_best(pop)

	# para individuo da população
	for i in range(len(pop)):

		# sorteia a probabilidade
		prob = np.random.randint(1,101)

		# testa a probabilidade
		if(prob< CONST_MUTATION_RATE):

			# escolhe os limites
			c1,c2 = np.random.choice(range(0,CONST_CITIES),2, replace=False)

			# salva a sublista em uma variavel aleatória
			aux = pop[i].cro[c1:c2]

			# inverte a sublista
			aux.reverse()

			# substitu a sublista pela invertida
			pop[i].cro[c1:c2] = aux

	# substitui o pior pelo melhor da população anterior
	pop[pop.index(get_worst(pop))] = copy.deepcopy(best)

	# retorn pop
	return pop

# Pega o melhor individuo
def get_best(pop):
	return min(pop,key = lambda x:x.fit)

# Pega o pior individuo
def get_worst(pop):
	return max(pop,key = lambda x:x.fit)

# plot as cidades
def plot_graph(city, m):

	# incializa variaveis com vazio
	x = []
	x1 = []
	y1 = []
	y = []

	# para cada cidade
	for i in range(len(CITIES)):
		# adicona os valoeres nas respectivas variaveis
		x.append(city[i][0])
		y.append(city[i][1])
		x1.append(CITIES[i][0])
		y1.append(CITIES[i][1])

	# adicona novamente o primeiro ponto para fechar o circuito
	x.append(city[0][0])
	y.append(city[0][1])
	x1.append(CITIES[0][0])
	y1.append(CITIES[0][1])

	# cria uma figura para salvar
	plt.figure()
	
	# plota a cidade passada como parametro
	plt.plot(x,y,'x-')

	# plota a ordem original com uma espessura de linha menor e outra cor
	plt.plot(x1,y1,'x--', linewidth=0.5)

	# plot um grid
	plt.grid(color='grey', linestyle='-', linewidth=0.1)

	# salva com o nome especifico dependendo dos parametros
	nameFigure = 'figure2_'+str(CONST_CROSS_RATE)+'_'+str(CONST_MUTATION_RATE) +'_'+str(CONST_POPULATION_SIZE)+'_'+str(CONS_GENERATIONS)+'_'+str(m)+'.png'

	# salva no disco
	plt.savefig(nameFigure)

	# mostra se necessário
	# plt.show()

# le arquivo do tipo tsp
def open_file(file):
	# Open input file
	
	# abre o arquivo em modo de leitura
	infile = open(file, 'r')

	# le as informações de cabecalho
	name = infile.readline().strip().split()[1] # NAME
	fileType = infile.readline().strip().split()[1] # TYPE
	comment = infile.readline().strip().split()[1] # COMMENT
	dimension = infile.readline().strip().split()[2] # DIMENSION
	edgeWeightType = infile.readline().strip().split()[1] # EDGE_WEIGHT_TYPE
	infile.readline()

	# lista de nós
	nodelist = []

	# para cidade presente no arquivo
	for i in range(int(dimension)):

			#  le os valores de coordena das cidades
			x,y = infile.readline().strip().split()[1:]

			# adiciona no formato desejado
			nodelist.append([int(x), int(y)])

	# fecha o arquivo
	infile.close()

	# atribui globalmente os valores de quantidade de cidades e vetor de cidades
	global	CONST_CITIES 
	CONST_CITIES = int(dimension)
	global	CITIES
	CITIES = nodelist

'''
MAIN
'''

# chama funcao para ler o arquivo
open_file('att48.tsp')


# nome de arquivo de log
nameFile = 'file2_'+str(CONST_CROSS_RATE)+'_'+str(CONST_MUTATION_RATE) +'_'+str(CONST_POPULATION_SIZE)+'_'+str(CONS_GENERATIONS)+'.txt'

# abre um aquivo no modo de escrita
f = open(nameFile, 'w+')

# para um determinado numero de iteracoes
for m in range(0,5):

	# escreve qual e a iterecao
	f.write('Exe: ' + str(m) + '\n')

	# gera a populacao inicial
	pop = [Cromossomo() for i in range(CONST_POPULATION_SIZE)]
	# # Mostra a populacao inicial
	print('--------------------------------')
	print('-------Populacao Inicial--------')
	show_population(pop)

	# inicia o melhor com vazio
	best = []

	# para todas as geracoes
	for i in range(0,CONS_GENERATIONS):		

		# realiza o torneio
		pop = make_tournament(pop)

		# mostra a populacao se necessario
		# show_population(pop)

		# realiza o obx
		# pop = make_obx(pop)
		pop = make_ox(pop)

		# mostra populacao se necessario
		# show_population(pop)

		# realiza mutacao por troc de elemento
		pop = make_ord_element_mutation(pop)
		# mostra populacao se necessario
		# show_population(pop)

		# caso queira selecoiona a mutacao por sublista
		# pop = make_inversion_sublist_mutation(pop)

		# mostra a populacao se necessario

		# show_population(pop)

		# A cada 100 geracoes mostra o melhor cromossomo
		if(i % 100 == 0 ):

			# salva o melhor daquela iterecao
			best = min(pop,key = lambda x:x.fit)

			# escreve o melhor o melhor no arquivo
			f.write('Generation: ' + str(i) + ' Best_Fit: '+str(best.fit)+'\n')

			# escrve o melhor na tela
			print('Generation: ' + str(i) + ' Best_Fit: '+str(best.fit))

	# plota o grafico para a primeira execucao
	plot_graph(best.cro,m)

	# Mostra a populacao final
	print('--------------------------------')
	print('-------Populacao Final----------')
	show_population(pop)

	# escreve no arquivo o melhor cromossomo
	f.write('\n\n\nBest: ')

	for k in range(len(best.cro)):
		f.write(str(best.cro[k]))
	f.write('\n\n\n')

# fecha o arquivo
f.close()