
"""
Using a connected dominating set (CDS) to serve as the virtual backbone of a wireless network
is an effective way to save energy and alleviate broadcasting storm. Since nodes may fail due to accidental damage
or energy depletion, it is desirable that the virtual backbone is fault tolerant. A node set C is an m-fold connected
dominating set (m-fold CDS) of graph G if every node in V(G)\C has at least m-neighbors in C and the subgraph of G
induced by C is connected. A greedy algorithm is presented to compute an m-fold CDS in a general graph, which has
size at most 2+ln(Delta + m - 2) times that of a minimum m-fold CDS, where Delta is the maximum degree of the graph.
This result improves on the previous best known performance ratio of 2*H(Delta + m - 1) for this problem, where H is
the harmonic number.

"""




""" 
A Python Class
A simple Python graph class, demonstrating the essential 
facts and functionalities of graphs.
"""
class Graph(object):

	def __init__(self, graph_dict={}):
		""" 
		Initializes a graph object
		"""
		self.__graph_dict = graph_dict

	def vertices(self):
		""" 
		Returns the vertices of a graph 
		"""
		return list(self.__graph_dict.keys())

	def edges(self):
		""" 
		Returns the edges of a graph 
		"""
		return self.__generate_edges()

	def add_vertex(self, vertex):
		""" 
		If the vertex "vertex" is not in 
		self.__graph_dict, a key "vertex" with an empty
		list as a value is added to the dictionary. 
		Otherwise nothing has to be done. 
		"""
		if vertex not in self.__graph_dict:
			self.__graph_dict[vertex] = []

	def add_edge(self, edge):
		""" 
		Assumes that edge is of type set, tuple or list; 
		between two vertices can be multiple edges! 
		"""
		edge = set(edge)
		vertex1 = edge.pop()
		if edge:
			# not a loop
			vertex2 = edge.pop()
		else:
			# a loop
			vertex2 = vertex1
		if vertex1 in self.__graph_dict:
			self.__graph_dict[vertex1].append(vertex2)
		else:
			self.__graph_dict[vertex1] = [vertex2]

	def __generate_edges(self):
		""" 
		A static method generating the edges of the 
		graph "graph". Edges are represented as sets 
		with one (a loop back to the vertex) or two 
		vertices 
		"""
		edges = []
		for vertex in self.__graph_dict:
			for neighbour in self.__graph_dict[vertex]:
				if {neighbour, vertex} not in edges:
					edges.append({vertex, neighbour})
		return edges

	def __str__(self):
		res = "vertices: "
		for k in self.__graph_dict:
			res += str(k) + " "
		res += "\nedges: "
		for edge in self.__generate_edges():
			res += str(edge) + " "
		return res

	def find_isolated_vertices(self):
		""" 
		Returns a list of isolated vertices. 
		"""
		graph = self.__graph_dict
		isolated = []
		for vertex in graph:
			#print(isolated, vertex)
			if not graph[vertex]:
				isolated += [vertex]
		return isolated


	def find_num_of_components(self):
		"""
		Return the number of componenets in a graph.. 
		"""
		graph = self.__graph_dict
		visited = []
		for vertex1 in graph:
			flag = 0
			for i in range(0, len(visited)):
				if vertex1 in visited[i] :
					flag = 1
			if flag:
				continue

			temp = []
			temp += vertex1
			for vertex2 in graph :
				#if vertex2 in visited:
				#	continue
				#if vertex1 == vertex2:
				#	continue
				if ( self.find_path(vertex1, vertex2) ):
					temp += vertex2

			visited.append(temp)

		return len(visited)


	def find_N(self, C):
		""" 
		For a node u e V , denote by N_c(u) the set of nodes in C which are adjacent with u in G 
		"""
		graph = self.__graph_dict
		N_c = {}
		for vertex1 in graph:
			#if vertex1 in C:
			#		continue
			all_paths = []
			temp = []
			for vertex2 in C:
				all_paths =  self.find_all_paths(vertex1, vertex2)

				for i in all_paths:
					if len(i) == 2:
						 temp +=  i[1]

			N_c[vertex1] = temp

		return N_c


	def v_sub_c(self, C):
		""" 
		V\C.. The elements in V that don't belong in C
		"""
		V =  self.__graph_dict
		subgraph = []

		for vertex1 in V:
			if vertex1 not in C:
				subgraph += vertex1

		return subgraph


	
	
	def V_i(self, C, m):
		""" 
		Returns: for i = 0...m-1, V_i = {u e V\C : |Nc| = i}.
		We use a dictionary of string numbers to enumerate the element inside the dictionary
		"""
		graph = self.__graph_dict

		V_sub_C = self.v_sub_c(C)

		V_i = {}

		N =  self.find_N(C)

		for i in range (0,m):
			V_i[str(i)] = []

		for vertex in V_sub_C:
			if  len(N[vertex]) > m -1:
				continue	
			V_i[str( len(N[vertex]) )] +=  vertex
		
		return V_i



	def V_m(self, C, m):
		""" 
		Returns: V_m = {u e V\C : |Nc| >= m}.
		No dictionaries and other shit here
		"""
		graph = self.__graph_dict

		V_sub_C = self.v_sub_c(C)

		V_m = []

		N =  self.find_N(C)

		for vertex in V_sub_C:
			if  len(N[vertex]) >= m: 
				V_m.append(vertex)
		return V_m


	def m(self, C, m):
		"""
		For all nodes in V, return a dictionary of type: ["u" :  m_c(u)].
		Where m_c(u) = 0, if u e C
			  m_c(u) = m-1, if u e V_0_C
			  m_c(u) = m-i, if u e V_i_C  (for i = 1,....,m)
		"""
		graph = self.__graph_dict
		m_c = {}

		vi = self.V_i(C, m)	
		vm = self.V_m(C, m)



		for vertex in graph:
			
			for el in vi:
				if int(el) >0  and (vertex in vi[el]):
					m_c[vertex] = m - int(el)

			if (vertex in C) or (vertex in vm):
				m_c[vertex] = 0
			if vertex in vi[str(0)]:
				m_c[vertex] = m-1
			
		return m_c


	def sum_m_C(self, dict):
		"""
		sum_m_C return the sum of m_c(u) for every u in V
		"""

		graph = self.__graph_dict
	
		sum = 0

		for vertex in graph:
			sum += dict[vertex]

		return sum



	def black_nodes(self, C):
		"""
		A node u is denoted black if u in C 
		"""
		graph = self.__graph_dict
		
		black_nodes = []

		for vertex in graph:
			if vertex in C:
				black_nodes.append(vertex)

		return black_nodes


	def gray_nodes(self, C, m):
		"""
		A node u is denoted gray if u in V_m_C
		"""
		graph = self.__graph_dict

		gray_nodes =  self.V_m(C, m)

		return gray_nodes 



	def red_nodes(self, C, m):
		"""
		A node u is denoted red if u in V_i_C (for all i= 1,....,m-1)
		"""

		graph = self.__graph_dict

		red_nodes = []  

		vi = self.V_i(C, m)

		for i in range(1, m):
			for j in vi[str(i)]:
				red_nodes.append(j)

		return red_nodes



	def white_nodes(self, C, m):
		"""
		A node u is denoted white if u in V_0_C
		"""

		graph = self.__graph_dict

		white_nodes =  []

		vi = self.V_i(C, m)	

		white_nodes = vi[str(0)]

		return white_nodes


	def spanning_subgraph(self, C):
		"""
		The spannign subgraph of G induced by 
		the edge set {e in E: e has at least on end in C}
		"""
		graph = self.__graph_dict

		spanning_sub = {}

		for vertex in graph:
			spanning_sub[vertex] = []

			if vertex in C:
				spanning_sub[vertex] = []
				spanning_sub[vertex] = g[vertex]
				continue
 
			flag =1
			for neighbour in g[vertex]:
				if (neighbour in C) and (flag ==1):
					spanning_sub[vertex] = []
					flag =0

				if neighbour in C:
					spanning_sub[vertex] += neighbour

		return spanning_sub


	def G_C(self, C):
		"""
		The subgraph of G induced by C
		"""
		graph = self.__graph_dict

		subgraph = {}

		for vertex in graph:
			#subgraph[vertex] = []  #remove comment to get all the 
			if vertex in C:
				subgraph[vertex] = []
				for neighbour in g[vertex]:
					if neighbour in C:
						subgraph[vertex] += neighbour
				
	
		return subgraph 


		
	def f_C(self, C, m):
		"""
		f_C = p(C) + q(C) + m(C) 
		"""		
		graph = self.__graph_dict		

		f_c = 0
		f_c = Graph(self.G_C(C)).find_num_of_components() + Graph(self.spanning_subgraph(C)).find_num_of_components() + self.sum_m_C(self.m(C, m))

		return f_c
		


	def delta_x_f_C(self, C, m, x):
		"""
		delta_x_f_C = f(C U {x}) - f(C)
		"""
		graph = self.__graph_dict

		delta = self.f_C( C + [x], m) - self.f_C( C, m)

		return delta




	def find_path(self, start_vertex, end_vertex, path=[]):
		""" 
		Find a path from start_vertex to end_vertex 
		in graph 
		"""
		graph = self.__graph_dict
		path = path + [start_vertex]
		if start_vertex == end_vertex:
			return path
		if start_vertex not in graph:
			return None
		for vertex in graph[start_vertex]:
			if vertex not in path:
				extended_path = self.find_path(vertex, 
											   end_vertex, 
											   path)
				if extended_path: 
					return extended_path
		return None
	

	def find_all_paths(self, start_vertex, end_vertex, path=[]):
		""" 
		Find all paths from start_vertex to 
		end_vertex in graph 
		"""
		graph = self.__graph_dict 
		path = path + [start_vertex]
		if start_vertex == end_vertex:
			return [path]
		if start_vertex not in graph:
			return []
		paths = []
		for vertex in graph[start_vertex]:
			if vertex not in path:
				extended_paths = self.find_all_paths(vertex, 
													 end_vertex, 
													 path)
				for p in extended_paths: 
					paths.append(p)
		return paths

	def is_connected(self, 
					 vertices_encountered = None, 
					 start_vertex=None):
		""" 
		Determines if the graph is connected 
		"""
		if vertices_encountered is None:
			vertices_encountered = set()
		gdict = self.__graph_dict        
		vertices = gdict.keys() 
		if not start_vertex:
			# chosse a vertex from graph as a starting point
			start_vertex = vertices[0]
		vertices_encountered.add(start_vertex)
		if len(vertices_encountered) != len(vertices):
			for vertex in gdict[start_vertex]:
				if vertex not in vertices_encountered:
					if self.is_connected(vertices_encountered, vertex):
						return True
		else:
			return True
		return False

	def vertex_degree(self, vertex):
		""" 
		The degree of a vertex is the number of edges connecting
		it, i.e. the number of adjacent vertices. Loops are counted 
		double, i.e. every occurence of vertex in the list 
		of adjacent vertices. 
		""" 
		adj_vertices =  self.__graph_dict[vertex]
		degree = len(adj_vertices) + adj_vertices.count(vertex)
		return degree

	def degree_sequence(self):
		""" 
		Calculates the degree sequence 
		"""
		seq = []
		for vertex in self.__graph_dict:
			seq.append(self.vertex_degree(vertex))
		seq.sort(reverse=True)
		return tuple(seq)

	@staticmethod
	def is_degree_sequence(sequence):
		""" Method returns True, if the sequence "sequence" is a 
			degree sequence, i.e. a non-increasing sequence. 
			Otherwise False is returned.
		"""
		# check if the sequence sequence is non-increasing:
		return all( x>=y for x, y in zip(sequence, sequence[1:]))
  

	def delta(self):
		""" the minimum degree of the vertices """
		min = 100000000
		for vertex in self.__graph_dict:
			vertex_degree = self.vertex_degree(vertex)
			if vertex_degree < min:
				min = vertex_degree
		return min
		
	def Delta(self):
		""" The maximum degree of the vertices """
		max = 0
		for vertex in self.__graph_dict:
			vertex_degree = self.vertex_degree(vertex)
			if vertex_degree > max:
				max = vertex_degree
		return max

	def density(self):
		""" method to calculate the density of a graph """
		g = self.__graph_dict
		V = len(g.keys())
		E = len(self.edges())
		return 2.0 * E / (V *(V - 1))

	def diameter(self):
		""" calculates the diameter of the graph """
		
		v = self.vertices() 
		pairs = [ (v[i],v[j]) for i in range(len(v)-1) for j in range(i+1, len(v))]
		smallest_paths = []
		for (s,e) in pairs:
			paths = self.find_all_paths(s,e)
			smallest = sorted(paths, key=len)[0]
			smallest_paths.append(smallest)

		smallest_paths.sort(key=len)

		# longest path is at the end of list, 
		# i.e. diameter corresponds to the length of this path
		diameter = len(smallest_paths[-1])
		return diameter

	@staticmethod
	def erdoes_gallai(dsequence):
		""" Checks if the condition of the Erdoes-Gallai inequality 
			is fullfilled 
		"""
		if sum(dsequence) % 2:
			# sum of sequence is odd
			return False
		if Graph.is_degree_sequence(dsequence):
			for k in range(1,len(dsequence) + 1):
				left = sum(dsequence[:k])
				right =  k * (k-1) + sum([min(x,k) for x in dsequence[k:]])
				if left > right:
					return False
		else:
			# sequence is increasing
			return False
		return True

   


if __name__ == "__main__":
	
	"""
	#C = []


	# testing C
	C= ["b", "e"]
	m = 2

	g = { 	"a" : ["b", "f", "d"],
			"b" : ["d", "e", "a"],
			"c" : ["d", "e", "f"],
			"d" : ["b", "c", "e", "a"],
			"e" : ["b", "c", "d"],
			"f" : ["a", "c"]
		}

	graph = Graph(g)


	# works like a charm
	print("N_C:")
	print ( graph.find_N(C) )

	print "V\C"
	print graph.v_sub_c( C)

	print "V_i"
	print graph.V_i(C, m)

	print "V_m"
	print graph.V_m(C, m)

	print "m========="
	print graph.m(C, m)


	print "SUM"
	print graph.sum_m_C(graph.m(C, m))	

	print "Black nodes"
	print graph.black_nodes(C)

	print "Gray nodes"
	print graph.gray_nodes(C, m)

	print "Red nodes"
	print graph.red_nodes(C, m)

	print "White nodes"
	print graph.white_nodes(C, m)

	print "Spanning subgraph"
	print graph.spanning_subgraph(C)

	print "Num of components in graph: p(c)"
	print Graph(graph.G_C(C)).find_num_of_components()

	print "Num of components in subgraph: q(c)"
	print Graph(graph.spanning_subgraph(C)).find_num_of_components()

	print "G_C"
	print graph.G_C(C)

	print "Noumero"
	print graph.f_C( C,m)

	print "delta_x_f_C"
	print graph.delta_x_f_C( C, m, "f")

	"""


	"""
	The algo
	"""

	"""
	test case 1 - [checked]
	g = { 	"a" : ["b", "f"],
			"b" : ["d", "e", "a"],
			"c" : ["d", "e"],
			"d" : ["b", "c", "e"],
			"e" : ["b", "c", "d"],
			"f" : ["a"]
		}

	m = 1
	C = [a,b,e]

	m = 2
	C = [a, b, c, e, f]
	"""

	"""
	test case 2 - [checked]

	g = { 	"a" : ["b", "c"],
			"b" : ["a", "c"],
			"c" : ["a", "b","d"],
			"d" : ["c", "e", "f", "g"],
			"e" : ["d", "f", "g"],
			"f" : ["e", "d", "g"],
			"g" : ["e", "d", "f", "k"],
			"k" : ["l", "g", "m","n", "h","i","j","p"],
			"h" : ["k", "j"],
			"i" : ["k", "j"],
			"j" : ["i", "h", "k", "o", "p"],
			"l" : ["k", "m"],
			"m" : ["l", "k", "n"],
			"n" : ["m", "k"],
			"o" : ["j", "p"],
			"p" : ["o", "j", "k", "q"],
			"q" : ["p", "s", "r", "t"],
			"s" : ["q", "r"],
			"r" : ["q", "s", "t"],
			"t" : ["r", "q"]

		}

	m = 1
	C =['c' , 'd' , 'g' , 'k' , 'q' , 'p']


	m = 2
	C = ['a', 'c', 'd', 'g', 'k', 'j', 'm', 'q', 'p', 'r']

	m = 3
	C = ['k', 'd', 'j', 'q', 'a', 'e', 'm', 'r', 'c', 'g', 'o', 'b', 'i', 'h', 'l', 'n', 'p', 's', 't']

	"""
	
	"""
	test case 3 - [checked]
	g = { 	"a" : ["b", "c", "e"],
			"b" : ["a", "e", "c","d"],
			"c" : ["a", "b","d","g"],
			"d" : ["c", "b", "f", "g"],
			"e" : ["a", "b", "f", "j"],
			"f" : ["e", "d", "g", "i","j"],
			"g" : ["c", "d", "f", "i", "h"],
			"i" : ["f", "g", "h","j"],
			"h" : ["g", "i","j"],
			"j" : ["h", "i","f","e"]

		}

	m = 1
	C =['e', 'f' ,'g']


	m = 2
	C = ['b', 'e', 'g', 'c', 'j']

	m = 3
	C = ['a','b', 'g', 'c', 'j','f','i' ]

	"""

	g = { 	"a" : ["b", "c"],
			"b" : ["a", "c"],
			"c" : ["a", "b","d"],
			"d" : ["c", "e", "f", "g"],
			"e" : ["d", "f", "g"],
			"f" : ["e", "d", "g"],
			"g" : ["e", "d", "f", "k"],
			"k" : ["l", "g", "m","n", "h","i","j","p"],
			"h" : ["k", "j"],
			"i" : ["k", "j"],
			"j" : ["i", "h", "k", "o", "p"],
			"l" : ["k", "m"],
			"m" : ["l", "k", "n"],
			"n" : ["m", "k"],
			"o" : ["j", "p"],
			"p" : ["o", "j", "k", "q"],
			"q" : ["p", "s", "r", "t"],
			"s" : ["q", "r"],
			"r" : ["q", "s", "t"],
			"t" : ["r", "q"]

		}

	graph = Graph(g)


	print "ALGORITHM"
	C = []
	m = 2

	for vertex in graph.v_sub_c( C):
		#print graph.delta_x_f_C( C, m, vertex)
		if -graph.delta_x_f_C( C, m, vertex) < 0:
			continue

		max_delta = 0
		for vertex2 in graph.v_sub_c( C):
		#	print vertex2,-graph.delta_x_f_C( C, m, vertex2)
			if (-graph.delta_x_f_C( C, m, vertex2)) > max_delta:
				max_delta = - graph.delta_x_f_C( C, m, vertex2)
				max_vertex = vertex2
				#print max_vertex

		if max_vertex not in C:
			C.append(max_vertex) 


	print C
