import networkx as nx
import random

ROOT_NAME = 'root'

TREE_FLAT = 1
TREE_BIFURCATING = 2
TREE_CASCADING = 3
TREE_RANDOM = 4

TOPOLOGY = (TREE_FLAT, TREE_BIFURCATING, TREE_CASCADING, TREE_RANDOM)


def get_topology_and_sizes(s):
    
    s = s.split()
    topology_name =  s[0]
    
    if topology_name == 'flat':
        topology = TREE_FLAT 
    elif topology_name == 'bifurcating':
        topology = TREE_BIFURCATING         
    elif topology_name  == 'cascading':
        topology = TREE_CASCADING         
    elif topology_name  == 'random':
        topology = TREE_RANDOM
    else:
        raise Exception(f'No such tree topology: "{topology_name}".')  
        
    sizes = list(map(int, s[1:]))
    
    if topology in (TREE_FLAT, TREE_CASCADING, TREE_RANDOM) and len(sizes) != 1:
        raise Exception(f'Only one number should be given to denote the size of the "{topology_name}" tree.')
        
    return {'topology': topology, 'sizes': sizes} 

def generate_tree(topology,  sizes):
    if topology in (TREE_FLAT, TREE_BIFURCATING):
        tree = generate_bifurcating_tree(sizes)
    elif topology == TREE_CASCADING:
        tree = generate_cascading_tree(sizes[0])
    elif topology == TREE_RANDOM:
        tree = generate_random(sizes[0])
    else:
        raise Exception(f'{topology} No such tree topology')
        
    return tree 

def generate_bifurcating_tree(sizes):
   
    G = nx.Graph() 
    G.add_node(ROOT_NAME)
    nodes = [ROOT_NAME]
    new_node = 0
    
    for size in sizes:
        
        new_nodes = []
        
        for node in nodes:
            
            for _ in range(size):
                G.add_edge(new_node, node)
                new_nodes.append(new_node)
                new_node += 1
                
        nodes = new_nodes
        
    return G

def generate_cascading_tree(size):
    
    G = nx.Graph()
    G.add_node('root')
    
    parent_node = 'root'
    new_node = 0
    
    for _ in range(size):
      
        G.add_edge(new_node, parent_node)
        new_node += 1
        G.add_edge(new_node, parent_node)
        
        parent_node = new_node
        new_node += 1
        
    return G

def generate_random_tree(size):
    
    G = nx.Graph()
    G.add_node(ROOT_NAME)
    nodes = [ROOT_NAME]
    new_node = 0
    
    for _ in range(size):
        
        random_parent = nodes[random.randrange(len(nodes))]
        
        G.add_edge(new_node, random_parent)
        nodes.append(new_node)
        
        new_node += 1
         
    return G

def save_tree(G, tree_fn):
    nx.write_gexf(G, tree_fn)

def read_tree(tree_fn):
    return nx.read_gexf(tree_fn)
