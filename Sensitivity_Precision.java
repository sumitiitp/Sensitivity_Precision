import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Arrays;


public class Sensitivity_Precision {

    public static void main(String[] args) {
        int n = 5;
        Naive_Approach na = new Naive_Approach();
        Hasse_Diagram_Approach hda = new Hasse_Diagram_Approach();
        
        long startTime = System.nanoTime();
        na.calculateSensitivity(n);
        long stopTime = System.nanoTime();
        long timeNaive = stopTime - startTime;
        System.out.println("Time taken by Naive Approach = " + timeNaive);
        
        startTime = System.nanoTime();
        hda.generateHasseDiagram(n);
        stopTime = System.nanoTime();
        long timeHasse = stopTime - startTime;
        System.out.println("Time taken by Hasse Diagram based Approach = " + timeHasse);
    }
}


class Precision {
    private int numerator;
    private int denominator;

    public Precision() {

    }
    
    public Precision(int numerator, int denominator) {
        this.numerator = numerator;
        this.denominator = denominator;
    }

    public Precision(Precision precision) {
        this.numerator = precision.numerator;
        this.denominator = precision.denominator;
    }
    
    public int getNumerator() {
        return numerator;
    }

    public int getDenominator() {
        return denominator;
    }

    public void setNumerator(int numerator) {
        this.numerator = numerator;
    }

    public void setDenominator(int denominator) {
        this.denominator = denominator;
    }

    @Override
    public String toString() {
        return "Precision{" + "numerator=" + numerator + ", denominator=" + denominator + '}';
    }
}


class SearchNode {
    private boolean isPresent;
    private int id;

    public SearchNode() {
    }

    public SearchNode(boolean isPresent, int id) {
        this.isPresent = isPresent;
        this.id = id;
    }

    public boolean isIsPresent() {
        return isPresent;
    }

    public int getId() {
        return id;
    }

    public void setIsPresent(boolean isPresent) {
        this.isPresent = isPresent;
    }

    public void setId(int id) {
        this.id = id;
    }

    @Override
    public String toString() {
        return "SearchNode{" + "isPresent=" + isPresent + ", id=" + id + '}';
    }
}


class TrieNode {
    int key; 
    TrieNode[] children; 
    boolean isLeaf;
    int Id;
  
    public TrieNode() { 
        this.isLeaf = false;
        this.Id = -1;
    } 

    public TrieNode(int key) { 
        this.key = key; 
        this.isLeaf = false;
        this.Id = -1;
    }   
    
    public TrieNode(int key, int size) { 
        this.key = key; 
        this.children = new TrieNode[size];
        this.isLeaf = false;
        this.Id = -1;
    } 

    public TrieNode(int key, int size, boolean isleaf, int id) { 
        this.key = key; 
        this.children = new TrieNode[size];
        this.isLeaf = true;
        this.Id = id;
    } 
    
    public TrieNode(TrieNode trieNode) { 
        this.key = trieNode.key; 
        this.children = trieNode.children.clone();
        this.isLeaf = trieNode.isLeaf;
        this.Id = trieNode.Id;
    } 
    
    @Override
    public String toString() {
        return "TrieNode{" + "key=" + key + ", children=" + children + ", isLeaf=" + isLeaf + ", Id=" + Id + '}';
    }
}


class Trie {
    TrieNode root;  // Root of Trie 
  
    // Constructor 
    Trie() {  
        this.root = null; 
    } 
    
    /* Insert partition P whose id is 'id' into the trie */
    void insert(int[] P, int n, int id) {
        int k = P.length; 
        int required_sum = n;
        int level = 0; 
        int index, max, min, noChildren;
        
        /* If the tree is empty, return a new node */
        if (this.root == null) { 
            int key = -1;
            max = required_sum - (k-level) + 1;
            min = (int)Math.ceil(required_sum/(k-level));
            noChildren = max - min + 1;
            this.root = new TrieNode(key, noChildren); // -1 to denote it is root
        } 
        TrieNode trieNode = root; 

        for (level = 0; level < k-1; level++)  {
            int key = P[level];
            min = (int)Math.ceil(required_sum/(k-level));
            index = key - min;
            
            required_sum = required_sum - key;
            if (trieNode.children[index] == null) {
                max = required_sum - (k-level-1) + 1;
                min = (int)Math.ceil(required_sum/(k-level-1));
                noChildren = max - min + 1;
                trieNode.children[index] = new TrieNode(key,noChildren);
            }
            trieNode = trieNode.children[index];   
        }
        int key = P[level];
        index = 0;
        trieNode.children[index] = new TrieNode(key,0, true, id);
    } 

    /* Search a partition P into the m-ary tree */
    SearchNode search(int[] P, int n) {
        int level; 
        int k = P.length; 
        int required_sum = n;
        int index; 
        
        TrieNode trieNode = root; 
       
        for (level = 0; level < k; level++) {
            //index = P[level]; 
            int key = P[level];
            int min = (int)Math.ceil(required_sum/(k-level));
            index = key - min;
            
            required_sum = required_sum - key;
            if (trieNode.children[index] == null) {
                SearchNode searchNode = new SearchNode(false, -1);
                return searchNode; 
            }
            trieNode = trieNode.children[index]; 
        } 
       
        if (trieNode == null) {
            SearchNode searchNode = new SearchNode(false, -1);
            return searchNode; 
        } else {
            SearchNode searchNode = new SearchNode(true, trieNode.Id);
            return searchNode; 
        }
    }  
}


class Vertex {
    private int vertexId;
    private int elements[];
    private int combinatorial;
    private boolean[] TP;
    
    public Vertex() {
        
    }

    public Vertex(int nodeId, int[] elements, int combinatorial, int size) {
        this.vertexId = nodeId;
        this.elements = elements.clone();
        this.combinatorial = combinatorial;
        this.TP = new boolean[size];
    }

    public Vertex(Vertex vertex) {
        this.vertexId = vertex.vertexId;
        this.elements = vertex.elements.clone();
        this.combinatorial = vertex.combinatorial;
        this.TP = vertex.TP.clone();
    }
    
    public int getVertexId() {
        return vertexId;
    }

    public int[] getElements() {
        return elements;
    }

    public int getCombinatorial() {
        return combinatorial;
    }

    public boolean[] getTP() {
        return TP;
    }

    public boolean getTP(int index) {
        return TP[index];
    }
    
    public void setVertexId(int vertexId) {
        this.vertexId = vertexId;
    }

    public void setElements(int[] elements) {
        this.elements = elements.clone();
    }

    public void setCombinatorial(int combinatorial) {
        this.combinatorial = combinatorial;
    }

    public void setTP(boolean[] TP) {
        this.TP = TP;
    }

    public void setTP(int index) {
        this.TP[index] = true;
    }
    
    @Override
    public String toString() {
        return "Vertex{" + "vertexId=" + vertexId + ", elements=" + Arrays.toString(elements) + ", combinatorial=" + combinatorial + ", TP=" + Arrays.toString(TP) + '}';
    }

}


class Naive_Approach {
    public List<List<List<Integer>>> obtainClusteringResults(int n) {
        int points[] = new int[n];
        for(int i = 0; i < n; i++) {
            points[i] = i+1;
        }
        List<List<List<Integer>>> results = new ArrayList<>();
        generateClusters(points, 0, new ArrayList<>(), results);    
        return results;
    }
    
     public static void generateClusters(int[] points, int index, List<List<Integer>> current, List<List<List<Integer>>> results) {
        if (index == points.length) {
            results.add(cloneList(current));
            return;
        }

        // Try to add the current point to each existing cluster
        for (int i = 0; i < current.size(); i++) {
            current.get(i).add(points[index]);
            generateClusters(points, index + 1, current, results);
            current.get(i).remove(current.get(i).size() - 1); // backtrack
        }

        // Create a new cluster with the current point
        List<Integer> newCluster = new ArrayList<>();
        newCluster.add(points[index]);
        current.add(newCluster);
        generateClusters(points, index + 1, current, results);
        current.remove(current.size() - 1); // backtrack
    }

    public static List<List<Integer>> cloneList(List<List<Integer>> list) {
        List<List<Integer>> clone = new ArrayList<>();
        for (List<Integer> item : list) {
            clone.add(new ArrayList<>(item));
        }
        return clone;
    }
    
	//Obtain the number of common elements in both the list. Elements in the list are arranged in ascending order.
    public int commonElement(List<Integer> list1, List<Integer> list2) {
        int count = 0;
        int index = 0;
        for(int i = 0; i < list1.size(); i++) {
            int a = list1.get(i);
            for(int j = index; j < list2.size(); j++) {
                if(a == list2.get(j)) {
                    count++;
                    index = j++;
                    break;
                }
            }
        }
        return count;
    }
    
    public int TP(int[][] contingencyTable, int row, int col) {
        int tp = 0;
        for(int i = 0; i < row; i++) {
            for(int j = 0; j < col; j++) {
                tp = tp + nC2(contingencyTable[i][j]);
            }
        }
        return tp;
    }
    
    public int nC2(int n) {
        if (n <= 1) {
            return 0;
        } else {
           return n*(n-1) / 2;
        }
    }
    
    public int FP(int[][] contingencyTable, int row, int col, int tp) {
        int fp = 0;
        for(int j = 0; j < col; j++) {
            int sum = 0;
            for(int i = 0; i < row; i++) {
                sum = sum + contingencyTable[i][j];
            }
            fp = fp + nC2(sum);
        }
        return fp - tp;
    }
    
    public void printMatrix(int[][] matrix, int row, int col) {
        for(int i = 0; i < row; i++) {
            for(int j = 0; j < col; j++) {
                System.out.print(matrix[i][j] + "   ");
            }
            System.out.println();
        }
    }
    
    public void calculateSensitivity(int n) {
        List<List<List<Integer>>> results = new ArrayList<>();
        results = obtainClusteringResults(n);
        
        /* System.out.println("The possible clustering results are...");
        for (List<List<Integer>> result : results) {
            System.out.println(result);
        } */
        
        int size = (n*(n-1))/2+1;
        boolean[][] matrix = new boolean[size][size];
        matrix[0][0] = true; /* Value of Precision: 0/0 */
        matrix[0][1] = true; /* Value of Precision: 0/1 */
        matrix[1][1] = true; /* Value of Precision: 1/1 */
        
        for (List<List<Integer>> result1 : results) {
            //System.out.println(result1);
            int row = result1.size();
            for (List<List<Integer>> result2 : results) {
                //System.out.println("    " + result2);
                int col = result2.size();
                int contingencyTable[][] = new int[row][col];
                for(int i = 0; i < row; i++) {
                    for(int j = 0; j < col; j++) {
                        contingencyTable[i][j] = commonElement(result1.get(i), result2.get(j));
                    }
                }
                /* System.out.println("Matrix is as follows...");
                printMatrix(contingencyTable, row, col); */
                
                int tp = TP(contingencyTable, row, col);
                int fp = FP(contingencyTable, row, col, tp);
                Precision precision = new Precision(tp, tp+fp);
                //System.out.println("precision = " + precision);
                
                if(precision.getNumerator() != 0) {
                    Precision reducedPrecision = new Precision();
                    reducedPrecision = reduceFraction(tp, tp+fp);
                    matrix[reducedPrecision.getNumerator()][reducedPrecision.getDenominator()] = true;
                }
            }
        }
        int sensitivity = 3; /* Value of Precision: 0/0, 1/1, 0/1 that is alreday considered. */
        /* Compute sensitivity */
        for(int i = 1; i < size; i++) {
            for(int j = i+1; j < size; j++) {
                if(matrix[i][j]) {
                    sensitivity++;
                }
            }
        }
        System.out.println("Sensitivity of Precision = " + sensitivity);
    }
    
    public Precision reduceFraction(int x, int y) {
        int d = gcd(x, y);  
        x = x / d;  
        y = y / d;  
        Precision pr = new Precision(x,y);
        return pr;
    }  
    
    public int gcd(int a, int b)  {
        if (b == 0) {  
            return a;  
        }
        return gcd(b, a % b);  
    }  
}


class Hasse_Diagram_Approach {

    public void generateHasseDiagram (int n) {
        // Calculate the number of partitions using Hardy-Ramanujan Asymptotic Partition Formula
        double pn_numerator = Math.PI * Math.sqrt((double)2*n/(double)3);
        double pn_denominator = 4*n* Math.sqrt(3);
        int pn = (int) Math.ceil(Math.exp(pn_numerator) / pn_denominator);
        
        int noVertices = 0;
        int noEdges = 0;
        
        /* Adjancy list representation of the graph (Considered the number of partitions given by 
         * Hardy-Ramanujan Asymptotic Partition Formula)
         * This will waste some space */ 
        LinkedList<Integer>[] adjList = new LinkedList[pn];
        /* An array of Vertex so that the partition can be accessed in constant time given the id 
         * for that corresponding partition */
        Vertex[] toAccessPartition = new Vertex[pn];
        for(int i = 0; i < pn; i++) {
            toAccessPartition[i] = new Vertex();
            adjList[i] = new LinkedList();
        }

        /* First node of the graph */
        int id = 0;
        int[] elements = {n};
        int combinatorial = findCombinatorial(elements);
        int size = (n*(n-1))/2+1;
        Vertex vertex = new Vertex(id, elements, combinatorial, size);
        
        toAccessPartition[id] = new Vertex(vertex);
        noVertices++;
        
        Queue<Integer> Q = new LinkedList();
        Q.add(id);

        Trie trie = new Trie();
        int lengthToChcekTrie = 0;
        while(!Q.isEmpty()) {
            Vertex curVertex = new Vertex();
            int currPartitionId = Q.poll();
            curVertex = toAccessPartition[currPartitionId];
            int len = curVertex.getElements().length;
            
            if(lengthToChcekTrie != len) {
                trie = new Trie();
            }
            if(len == n) { // All the partition have been generated
                break;
            }
            
            /* Generate all the partition of length len+1 from a partition curVertex.getElements() of 
             * length len */
            for(int i = 0; i < len; i++) {
                if ((len==1) || i == 0 || (curVertex.getElements()[i] != curVertex.getElements()[i-1]) ) {
                    int[] newPartition = new int[len + 1];
                    if(i-1 >= 0) {
                        System.arraycopy(curVertex.getElements(), 0, newPartition, 0, i-1+1);
                    }
                    System.arraycopy(curVertex.getElements(), i+1, newPartition, i, len-(i+1));
                    for(int j = 1; j <= curVertex.getElements()[i]/2; j++) {
                        newPartition[newPartition.length-2] = j;
                        newPartition[newPartition.length-1] = curVertex.getElements()[i]-j;
                        int[] sortedNewPartition = new int[newPartition.length];
                        sortedNewPartition = countingSort(newPartition);
                        
                        /* Chcek whether the generated partition has been alreday generated or not
                         using m-ary tree based data structure */
                        SearchNode searchNode = new SearchNode();
                        if(trie.root == null) {
                            searchNode = new SearchNode(false, -1);
                        } else {
                            searchNode = trie.search(sortedNewPartition, n);
                        }
                        noEdges++;
                        if(searchNode.isIsPresent()) {
                           int existingNodeid = searchNode.getId();
                           adjList[currPartitionId].add(existingNodeid);   
                        } else {
                            id++;   
                            combinatorial = findCombinatorial(sortedNewPartition);
                            vertex = new Vertex(id, sortedNewPartition, combinatorial, size);
                            adjList[currPartitionId].add(id);   
                            toAccessPartition[id] = new Vertex(vertex);
                            trie.insert(sortedNewPartition, n, id);
                            noVertices++;
                            Q.add(id);
                        }
                    }
                }
            }  
            lengthToChcekTrie = len;
        }
        /* System.out.println("Statistics about Hasse Diagram");
        System.out.println("    n = " + n);
        System.out.println("    No. of vertices = " + noVertices);
        System.out.println("    No. of edges = " + noEdges); */
        
        boolean[][] matrix = new boolean[size][size];
        matrix[0][0] = true; /* Value of Precision: 0/0 */
        matrix[0][1] = true; /* Value of Precision: 0/1 */
        matrix[1][1] = true; /* Value of Precision: 1/1 */
        for(int i = noVertices-2; i >= 0; i--) {
            combinatorial = toAccessPartition[i].getCombinatorial();
            toAccessPartition[i].setTP(combinatorial);
            for(int j = 0; j < adjList[i].size(); j++) {
                toAccessPartition[i].setTP(logicalOR(toAccessPartition[i].getTP(), toAccessPartition[adjList[i].get(j)].getTP()));
            }
            for(int j = 1; j < toAccessPartition[i].getTP().length; j++) {
                Precision pr = new Precision();
                if(toAccessPartition[i].getTP(j) == true) {
                    pr = reduceFraction(j, toAccessPartition[i].getCombinatorial());
                    matrix[pr.getNumerator()][pr.getDenominator()] = true;
                }
            }
        }
       
        
        int sensitivity = 3; /* Value of Precision: 0/0, 1/1, 0/1 that is alreday considered. */
        /* Compute sensitivity */
        for(int i = 1; i < size; i++) {
            for(int j = i+1; j < size; j++) {
                if(matrix[i][j]) {
                    sensitivity++;
                }
            }
        }
        System.out.println("Sensitivity of Precision = " + sensitivity);
    }

    public boolean[] logicalOR(boolean[] A, boolean[] B) {
        for(int i = 0; i < A.length; i++) {
            A[i] = A[i] | B[i];
        }
        return A;
    }
    
    public int findCombinatorial(int[] array) {
        int combinatorial = 0;
        for(int i = 0; i < array.length; i++) {
            if(array[i] >= 2) {
                combinatorial = combinatorial + (array[i]*(array[i]-1))/2;
            }
        }
        return combinatorial;
    }
    
    public Precision reduceFraction(int x, int y) {
        int d = gcd(x, y);  
        x = x / d;  
        y = y / d;  
        Precision pr = new Precision(x,y);
        return pr;
    }  
  
    
    public int gcd(int a, int b)  {
        if (b == 0) {  
            return a;  
        }
        return gcd(b, a % b);  
    } 
    
    public int[] countingSort(int[] array) { 
        int[] aux = new int[array.length];

        // Find the smallest and the largest value in the array
        int min = array[0];
        int max = array[0];
        for (int i = 1; i < array.length; i++) {
            if (array[i] < min) {
                min = array[i];
            } else if (array[i] > max) {
                max = array[i];
            }
        }

        // Initialize array of frequencies
        int[] counts = new int[max - min + 1];

        // Initialize the frequencies
        for (int i = 0;  i < array.length; i++) {
          counts[array[i] - min]++;
        }

        // Recalculate the array - create the array of occurences
        counts[0]--;
        for (int i = 1; i < counts.length; i++) {
          counts[i] = counts[i] + counts[i-1];
        }

        for (int i = array.length - 1; i >= 0; i--) {
            aux[counts[array[i] - min]--] = array[i];
        }
        
        int[] auxDescending = new int[array.length];
        for (int i = 0; i < array.length; i++) {
            auxDescending[i] = aux[array.length - i - 1];
        }
        return auxDescending;
    } 
}
