import numpy as np

def creer_matrice_tridiagonale_inferieure_mod_p(n, element_sous_diag, p):
    matrice = np.zeros((n, n), dtype=int)
    
    np.fill_diagonal(matrice, 1)
    
    for i in range(1, n):
        matrice[i, i-1] = element_sous_diag  
    
    matrice_mod_p = matrice % p

    return matrice_mod_p

def creer_matrice_tridiagonale_superieure_mod_p(n, element_super_diag, p):
    matrice = np.zeros((n, n), dtype=int)
    
    np.fill_diagonal(matrice, 1)
    
    for i in range(n-1):
        matrice[i, i+1] = element_super_diag  
    matrice_mod_p = matrice % p

    return matrice_mod_p

def calcul_determinant(matrice):
    if matrice.shape[0] != matrice.shape[1]:
        return "Error: The Matrix should be square."
    
    try:
        determinant = np.linalg.det(matrice)
        return determinant
    except np.linalg.LinAlgError:
        return "Error while calculating the determinant."

def test_determinant_un(matrice):
    determinant = calcul_determinant(matrice)
    
    if isinstance(determinant, (int, float)):
        return np.isclose(determinant, 1)
    else:
        return False

def matrix_power_mod_p(A, n, p):
    if A.shape[0] != A.shape[1]:
        raise ValueError("The Matrix should be square.")
    
    if n == 0:
        return np.eye(A.shape[0], dtype=int) % p
    
    if n == 1:
        return A % p
    
    result = np.eye(A.shape[0], dtype=int)  
    for _ in range(n):
        result = np.dot(result, A) % p  
    
    return result

def matrice_mod_p(matrice, p):
    
    matrice_mod = np.mod(matrice, p)
    
    return matrice_mod

def matrix_inverse(A):
    if A.shape[0] != A.shape[1]:
        raise ValueError("The Matrix should be square.")
    
    try:
        A_inv = np.linalg.inv(A)
        return A_inv
    except np.linalg.LinAlgError:
        return "The Matrix is not invertible."

def invp(A,p):
    return(matrice_mod_p(matrix_inverse(A), p))

def is_prime(num):
    
    if num <= 1:
        return False
    for i in range(2, int(num**0.5) + 1):
        if num % i == 0:
            return False
    return True

def a_b1(n):
    if n == 3:
        a = 4  
        b = 2 
        l=4
        return a, b ,l
    else:
        q = None
        for candidate in range(2, n):  
            if is_prime(candidate) and n % candidate == 1:
                q = candidate
                break
        if q is None:
            return "No prime number q has been found."
    k = 1
    l = q ** k + 1
    while l < 3 * (n - 1):
        k += 1
        l = q ** k + 1
        a = q + 1 
        b = a + q 

    return a,b,l

def create_application1(n,p):
    A=matrix_power_mod_p(creer_matrice_tridiagonale_superieure_mod_p(n,a_b1(n)[0],p), a_b1(n)[2],p)
    B=matrix_power_mod_p(creer_matrice_tridiagonale_inferieure_mod_p(n,a_b1(n)[1],p),a_b1(n)[2],p)
    C=invp(A,p)
    D=invp(B,p)
    

    app = {
        1: B,  
        2: C,  
        3: D,  
        
    }
    return app

def s1(app, index):
    if index in app:
        return(app[index]) 
    else:
        return(f"The index {index} doesn't exxist in the application.")

def create_application2(n,p):
    A=matrix_power_mod_p(creer_matrice_tridiagonale_superieure_mod_p(n,a_b1(n)[0],p), a_b1(n)[2],p)
    B=matrix_power_mod_p(creer_matrice_tridiagonale_inferieure_mod_p(n,a_b1(n)[1],p),a_b1(n)[2],p)
    C=invp(A,p)
    D=invp(B,p)

    app = {
        1: A,  
        2: C, 
        3: D,  
        
    }
    return app

def s2(app, index):
    if index in app:
        return(app[index])  
    else:
        return(f"The index {index} doesn't exist in the application.")

def create_application3(n,p):
    A=matrix_power_mod_p(creer_matrice_tridiagonale_superieure_mod_p(n,a_b1(n)[0],p), a_b1(n)[2],p)
    B=matrix_power_mod_p(creer_matrice_tridiagonale_inferieure_mod_p(n,a_b1(n)[1],p),a_b1(n)[2],p)
    C=invp(A,p)
    D=invp(B,p)
    
    app = {
        1: A,  
        2: D,  
        3: B,  
        
    }
    return app

def s3(app, index):
    if index in app:
        return(app[index])  
    else:
        return(f"The index {index} doesn't exist in the application.")

def create_application4(n,p):
    A=matrix_power_mod_p(creer_matrice_tridiagonale_superieure_mod_p(n,a_b1(n)[0],p), a_b1(n)[2],p)
    B=matrix_power_mod_p(creer_matrice_tridiagonale_inferieure_mod_p(n,a_b1(n)[1],p),a_b1(n)[2],p)
    C=invp(A,p)
    D=invp(B,p)
    
    app = {
        1: A,  
        2: C,  
        3: B,  
        
    }
    return app

def s4(app, index):
    if index in app:
        return(app[index])  
    else:
        return(f"The index {index} doesn't exist in the application.")

def trouver_matrice(M,n,p):
    A=matrix_power_mod_p(creer_matrice_tridiagonale_superieure_mod_p(n,a_b1(n)[0],p), a_b1(n)[2],p)
    B=matrix_power_mod_p(creer_matrice_tridiagonale_inferieure_mod_p(n,a_b1(n)[1],p),a_b1(n)[2],p)
    C=invp(A,p)
    D=invp(B,p)
    matrices = [A, B, C, D]
    
    for i, matrice in enumerate(matrices):
        if np.array_equal(invp(M,p), matrice):
            return i + 1  
    return 0  

def hash1(m, n,p):
    L = []
    E = s1(create_application1(n,p), m[0])
    L.append(E)
    
    for i in range(1, len(m)):  
        if trouver_matrice(E, n,p) == 1:
            E = s1(create_application1(n,p), m[i])
            L.append(E)
        elif trouver_matrice(E, n,p) == 2:
            E = s2(create_application2(n,p), m[i])
            L.append(E)
        elif trouver_matrice(E, n,p) == 3:
            E = s3(create_application3(n,p), m[i])
            L.append(E)
        elif trouver_matrice(E, n,p) == 4:
            E = s4(create_application4(n,p), m[i])
            L.append(E)

    return L

def produit_matrices(m,n,p):
    # initialization for the product!!
    produit = hash1(m, n,p)[0]
    
    for i in range(1, len( hash1(m, n,p))):
        produit = np.dot(produit,  hash1(m, n,p)[i])
    
    return matrice_mod_p(produit, p)

#example
produit_matrices([2,3,1,2], 4, 12157)
