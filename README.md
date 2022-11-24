# PHH
---
Persistent hypergraph homology algorithm
---

# How to use
```python
'''
Given a filtered hypergraph H = [ [[1],0], [[1,2],1], [[2,3],2], [[1,3],3]  ] where 
hyperedge [1] has value 0, [1,2] has value 1, [2,3] has value 2, [1,3] has value 3. 
'''
H = [ [ [1],0 ], [ [1,2],1 ], [ [2,3],2 ], [ [1,3],3 ]  ]
a = Hypergraph(H)
dim = 3  # homology dimension dim - 1
bar_inf = a.persistence_from_infimum(dim)   # get barcode by infimum chain complex
print(bar_inf)
bar_sup = a.persistence_from_supremum(dim)  # get barcode by supremum chain complex, note that inf and sup will give the same results.
print(bar_sup)

'''
Output should like this:
{'0': [[0, -1]], '1': [[3, -1]], '2': []}
{'0': [[0, -1]], '1': [[3, -1]], '2': []}
'''
```

