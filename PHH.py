# -*- coding: utf-8 -*-
"""
@author: liuxiang
persistent hypergraph homology
"""



class Hypergraph():
     
    @classmethod
    def get_face(cls,path,index):
        res = []
        for i in range(len(path)):
            if i!=index:
                res.append(path[i])
        return res
    
    @classmethod
    def get_max(cls,vector):
        if len(vector)==0:
            return -1
        else:
            return max(vector)
        
    @classmethod
    def add_two_column(cls,a,b):
        for item in a:
            if item in b:
                b.remove(item)
            else:
                b.add(item)
        return b
    
    @classmethod
    def get_copy_vector(cls,v):
        res = []
        for i in v:
            res.append(i)
        return set(res)
    
    @classmethod
    def solve_equation(cls,vector,B,L):
        if len(vector)==0:
            return [1,[]]
        solution = []
        N = Hypergraph.get_max(vector)
        while N!=-1 and L[N]!=-1:
            vector = Hypergraph.add_two_column(B[L[N]],vector)
            solution.append(L[N])
            N = Hypergraph.get_max(vector)
        if N==-1:
            return [1,solution]
        elif L[N]==-1:
            return [0,vector]
        
    @classmethod
    def combine_D_M(cls,D,D_L,M,length,MM):
        DM = []
        D_end = len(D)
        DM_L = []
        add_index = []
        for item in D:
            temp1 = []
            for i in item:
                temp1.append(i)
            DM.append(set(temp1))
        for i in range(length):
            DM_L.append(-1)
        for i in range(len(D_L)):
            DM_L[i] = D_L[i]
        
        c = 0
        zero_number = 0
        for j in range(len(M)):
            N = Hypergraph.get_max(M[j])
            while N!=-1 and DM_L[N]!=-1:
                M[j] = Hypergraph.add_two_column(DM[ DM_L[N] ], M[j])
                if DM_L[N]>=D_end:
                    MM[j] = Hypergraph.add_two_column(MM[ DM_L[N]-D_end+zero_number ], MM[j])
                N = Hypergraph.get_max(M[j])
                
            if N==-1:
                zero_number = zero_number + 1
                continue
            elif DM_L[N]==-1:
                DM_L[N] = c + D_end
                c = c + 1
                DM.append(M[j])
                add_index.append(j)
        return [ [DM,DM_L],D_end,add_index ]
    
    @classmethod
    def get_smaller_number(cls,a,ls):
        res = 0
        for i in ls:
            if i<a:
                res = res + 1
        return res
    
    
    def __init__(self,hyperedge):
        #self.Vertex = vertex # vertex
        #self.Vertex_Number = len(vertex)
        self.WH = hyperedge
        self.H = []
        self.KH = []
        self.W = []
        self.H_Number = 0
        self.KH_Number = 0
        self.set_KH()
        self.Sup_Basis = []
        self.Inf_Basis = []
        self.Sup_Matrix = []
        self.Inf_Matrix = []
        self.Sup_Weight = []
        self.Inf_Weight = []
        self.Sup_Dimension = []
        self.Inf_Dimension = []
        self.Sup_L = []
        self.Inf_L = []
        self.Inf_Index = []
        self.Diagram = {}
        
        
    def set_KH(self):
        for item in self.WH:
            self.H.append(item[0])
            self.W.append(item[1])
            self.KH.append(item[0])
        for item in self.H:
            for i in range(len(item)):
                face = Hypergraph.get_face(item,i)
                if face not in self.H:
                    if len(face)>0:
                        self.KH.append(face)
        self.H_Number = len(self.H)
        self.KH_Number = len(self.KH)
    
    
    def get_boundary(self,temp):
        if len(temp)==1:
            return set()
        res = set()
        for i in range(len(temp)):
            face = Hypergraph.get_face(temp,i)
            index = self.KH.index(face)
            res.add(index)
        return res
    
    
    def supremum_chain_complex(self):
        for i in range(self.KH_Number):
            self.Sup_L.append(-1)
        for j in range(self.H_Number):
            vector_now = set([self.KH.index(self.H[j])])
            solution1 = Hypergraph.solve_equation(vector_now,self.Sup_Basis,self.Sup_L)
            if solution1[0]==1:
                continue
            else:
                temp_boundary = self.get_boundary(self.H[j])
                solution2 = Hypergraph.solve_equation(temp_boundary,self.Sup_Basis,self.Sup_L)
                if solution2[0]==1:
                    self.Sup_Basis.append(solution1[1])
                    temp_N = len(self.Sup_Basis)
                    self.Sup_L[max(solution1[1])] = temp_N - 1
                    self.Sup_Matrix.append(set(solution2[1]))
                    self.Sup_Weight.append(self.W[j])
                    self.Sup_Dimension.append(len(self.H[j])-1)
                else:
                    self.Sup_Basis.append(solution2[1])
                    temp_N = len(self.Sup_Basis)
                    self.Sup_L[max(solution2[1])] = temp_N -1
                    self.Sup_Basis.append(solution1[1])
                    self.Sup_L[max(solution1[1])] = temp_N
                    self.Sup_Matrix.append(set())
                    temp_N = len(self.Sup_Matrix)
                    self.Sup_Matrix.append(set([temp_N-1]))
                    self.Sup_Weight.append(self.W[j])
                    self.Sup_Weight.append(self.W[j])
                    self.Sup_Dimension.append(len(self.H[j])-2)
                    self.Sup_Dimension.append(len(self.H[j])-1)
        
        
    def set_infimum_index(self):
        D = []
        partial_D = []
        D_L = []
        for i in range(self.KH_Number):
            D_L.append(-1)
        for j in range(self.H_Number):
            sigma_now = self.H[j]
            vector_now = set([self.KH.index(sigma_now)])
            temp_boundary = self.get_boundary(sigma_now)
            
            now_D = []
            now_partial_D = []
            now_D_L = []
            for item in D:
                now_D.append(Hypergraph.get_copy_vector(item))
            for item in partial_D:
                now_partial_D.append(Hypergraph.get_copy_vector(item))
            for item in D_L:
                now_D_L.append(item)
            [DM,D_end,partial_D_index] = Hypergraph.combine_D_M(now_D,now_D_L,now_partial_D,self.KH_Number,now_D)
            temp_boundary2 = Hypergraph.get_copy_vector(temp_boundary)
            # partial_sigma
            solution1 = Hypergraph.solve_equation(temp_boundary2,DM[0],DM[1])
            if solution1[0]==1:
                self.Inf_Index.append(j)
            # sigma
            sigma = Hypergraph.get_copy_vector(vector_now)
            solution2 = Hypergraph.solve_equation(sigma, DM[0], DM[1])
            if solution2[0]==1:
                temp_M = partial_D_index[max(solution2[1])-D_end]
                self.Inf_Index.append(temp_M)

            D.append(vector_now)
            temp_N = len(D)
            D_L[max(vector_now)] = temp_N - 1 
            partial_D.append(temp_boundary)

    
    def infimum_chain_complex(self):
        self.supremum_chain_complex()
        self.set_infimum_index()      
        candidate = []
        for i in range(len(self.Sup_Basis)):
            temp = max(self.Sup_Basis[i])
            if temp in self.Inf_Index:
                candidate.append(i)
        inf_candidate = []
        for i in range(len(self.Sup_Basis)):
            if i not in candidate:
                inf_candidate.append(i)
        for i in range(len(self.Sup_Matrix)):
            if i in candidate:
                temp1 = self.Sup_Matrix[i]
                temp2 = set()
                for item in temp1:
                    if item in candidate:
                        temp2.add(item)
                temp3 = set()
                for item in temp2:
                    temp4 = Hypergraph.get_smaller_number(item,inf_candidate)
                    temp5 = item - temp4
                    temp3.add(temp5)
                self.Inf_Matrix.append(temp3)
                self.Inf_Weight.append(self.Sup_Weight[i])
                self.Inf_Dimension.append(self.Sup_Dimension[i])
           
        
    def persistence_from_supremum(self,Dim):
        self.supremum_chain_complex()
        inf_bar = []
        finite_bar = []
        t = len(self.Sup_Matrix)
        L = []
        for i in range(t):
            L.append(-1)
        for j in range(t):
            N = Hypergraph.get_max(self.Sup_Matrix[j])
            while N!=-1 and L[N]!=-1:
                self.Sup_Matrix[j] = Hypergraph.add_two_column(self.Sup_Matrix[L[N]],self.Sup_Matrix[j])
                N = Hypergraph.get_max(self.Sup_Matrix[j])
            if N==-1:
                inf_bar.append(j)
            elif L[N]==-1:
                L[N] = j
        
        for j in range(t):
            low_j = Hypergraph.get_max(self.Sup_Matrix[j])
            if low_j!=-1:
                inf_bar.remove(low_j)
                finite_bar.append([low_j,j])
        
    
        for i in range(Dim):
            self.Diagram[str(i)] = []
        for i in inf_bar: 
            self.Diagram[str(self.Sup_Dimension[i])].append([self.Sup_Weight[i],-1])
        for item in finite_bar:
            if self.Sup_Weight[item[1]]>self.Sup_Weight[item[0]]:
                self.Diagram[str(self.Sup_Dimension[item[0]])].append([ self.Sup_Weight[item[0]],self.Sup_Weight[item[1]] ])
        return self.Diagram
    
    
    def persistence_from_infimum(self,Dim):
        self.infimum_chain_complex()
        inf_bar = []
        finite_bar = []
        t = len(self.Inf_Matrix)
        L = []
        for i in range(t):
            L.append(-1)
        for j in range(t):
            N = Hypergraph.get_max(self.Inf_Matrix[j])
            while N!=-1 and L[N]!=-1:
                self.Inf_Matrix[j] = Hypergraph.add_two_column(self.Inf_Matrix[L[N]],self.Inf_Matrix[j])
                N = Hypergraph.get_max(self.Inf_Matrix[j])
            if N==-1:
                inf_bar.append(j)
            elif L[N]==-1:
                L[N] = j
        
        for j in range(t):
            low_j = Hypergraph.get_max(self.Inf_Matrix[j])
            if low_j!=-1:
                inf_bar.remove(low_j)
                finite_bar.append([low_j,j])
        
    
        for i in range(Dim):
            self.Diagram[str(i)] = []
        for i in inf_bar: 
            self.Diagram[str(self.Inf_Dimension[i])].append([self.Inf_Weight[i],-1])
        for item in finite_bar:
            if self.Inf_Weight[item[1]]>self.Inf_Weight[item[0]]:
                self.Diagram[str(self.Inf_Dimension[item[0]])].append([ self.Inf_Weight[item[0]],self.Inf_Weight[item[1]] ])
        return self.Diagram
    
                    
            
    
    


    
    
    