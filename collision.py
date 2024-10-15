import numpy as np

def finding_collisions(num_inputs, n, p):
    hash_results = {}  
    collisions = [] 

    
    for i in range(num_inputs):
        m = tuple(np.random.randint(1, 4, size=n)) 
        hash_value = tuple(map(tuple, produit_matrices(m, n, p)))  

       
        if hash_value in hash_results:
            collisions.append((hash_results[hash_value], m)) 
        else:
            hash_results[hash_value] = m  
        
        if i % 100 == 0:  
            print(f"Processed {i} inputs...")

    
    if collisions:
        print(f"Collisions found: {len(collisions)}")
        for input_seq1, input_seq2 in collisions:
            print(f"Collision: Input sequence 1: {input_seq1}, Input sequence 2: {input_seq2}")
    else:
        print(f"No collisions found after checking {num_inputs} inputs.")

# Example
p = 12157  
n = 6   
num_inputs = 100  

np.random.seed(42)

finding_collisions(num_inputs, n, p)
