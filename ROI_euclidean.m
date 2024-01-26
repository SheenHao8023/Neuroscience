points = [
69.98 	-25.53 	-2.36 
64.29 	-57.11 	24.40 
55.64 	-8.89 	50.68 
51.78 	32.26 	31.21 
18.60 	45.84 	48.09 
23.13 	66.21 	12.41 
-20.35 	66.34 	14.41 
-16.96 	45.82 	49.08 
-50.09 	31.52 	35.45 
-51.70 	-10.50 	55.83 
-64.39 	-55.12 	27.81 
-70.03 	-25.74 	5.64 
];
distances = pdist(points, 'euclidean');
distance_matrix = squareform(distances);
writematrix(distance_matrix, 'distance_matrix.xlsx');