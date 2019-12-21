import java.io.FileWriter;
import java.io.PrintWriter;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.text.DecimalFormat;
import java.util.Arrays;

public class TravelingSalesman {
    static ThreadMXBean bean = ManagementFactory.getThreadMXBean( );

    /* define constants */
    static int numberOfTrials = 5;
    static int MAXINPUTSIZE = 13;
    static int MININPUTSIZE = 4;

    //set up variable to hold folder path and FileWriter/PrintWriter for printing results to a file
    static String ResultsFolderPath = "/home/alyssa/Results/TSP/"; // pathname to results folder 
    static FileWriter resultsFile;
    static PrintWriter resultsWriter;

    public static void main (String [] args)
    {
        double costMatrix [][] = GenerateRandomCircularGraphCostMatrix(10, 100);
        // print cost matrix that is used for test input
/*      System.out.println("\nCost Matrix:");
        System.out.println("============");
        for(int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++)
                System.out.printf("%11.6f   ", costMatrix[i][j]);
            System.out.println();
        }
        System.out.println();                                                                                               */

        // run the greedy TSP algorithm on the generated cost matrix
/*      Path greedyPath = GreedyTSP(costMatrix);
        System.out. print("\nGreedy Path   : ");
        for(int i = 0; i < 10; i++)
            System.out.printf("%2d   ", greedyPath.path[i]);
        System.out.printf("\nGreedy Cost   = %11.6f", greedyPath.pathCost);

        // run the brute force TSP algorithm on the generated cost matrix
        Path brutePath = bruteForceTSP(costMatrix);
        System.out. print("Brute Force Path : ");
        for(int i = 0; i < 10; i++)
            System.out.print(brutePath.path[i] + "  ");
        System.out.printf("\nBrute Force Cost = %11.6f", brutePath.pathCost);                                               */

        // run the ant colony TSP algorithm on the generated cost matrix
/*      Path antColPath = AntColonyTSP(costMatrix, (int)(costMatrix[0].length * 0.5),1, 150);
        System.out. print("\nAnt Col Force Path : ");
        for(int i = 0; i < 10; i++)
            System.out.printf("%2d   ", antColPath.path[i]);
        System.out.printf("\nAnt Col Cost =  %11.6f", antColPath.pathCost);                                                 */

        // run full experiment 3 times
        for( int run = 1; run <= 3; run++)
        {
            runFullExperiment("TSPBruteForce--Run" + run + ".txt");
           // runFullExperiment("TSPGreedy--Run" + run + ".txt");
           // runFullExperiment("TSPAntColony--Run" + run + ".txt");
           // runFullExperiment("HeuristicSolutions--Run" + run + ".txt");
        }
    }

    //brute force TSP wrapper function
    static Path bruteForceTSP(double [][] costMatrix)
    {
        int numVertices = costMatrix[0].length;

        int basePath [] = new int [numVertices+1];
        // set the base path
        for(int i = 0; i < numVertices; i++ )
            basePath[i] = i;
        basePath[numVertices] = 0;

        // create a new Path with the maximum value for pathCost
        Path best = new Path (numVertices);
        best.pathCost = Double.MAX_VALUE;

        // call the BruteForceTSP function with the path object to store the best path, costMatrix, the base path created
        bruteForceTSP(best, costMatrix, basePath, 1, numVertices);
        return best;
    }

    // brute force TSP function
    static void bruteForceTSP(Path best, double [][]costMatrix, int path[], int a, int b)
    {
        // once a and b indexes meet up
        if( a == b)
        {
            double cost = 0;
            // calculate the cost of the path created
            for(int i = 0; i < path.length-1; i++)
                cost += costMatrix[path[i]][path[i+1]];
            // if the cost is less than the current cost for the best path, update the best path
            if (cost < best.pathCost)
            {
                best.pathCost = cost;
                best.path = path;
            }
        }
        else
            // for each a index less than b index
            for(int i = a; i < b; i++)
            {
                // swap vertices at a and i
                path = swap(path, a, i);
                // recursively call bruteForceTSP
                bruteForceTSP(best, costMatrix, path, a+1, b);
                // swap vertices at a and i
                path = swap(path, a, i);
            }
    }

    // swap two values in an array
    static int [] swap(int [] path, int v1, int v2)
    {
        int [] newPath = path.clone();
        newPath[v1] = path[v2];
        newPath[v2] = path[v1];
        return newPath;
    }

    // greedy algorithm for TSP
    static Path GreedyTSP(double [][] costMatrix)
    {
        int numVertices = costMatrix[0].length;
        Path shortest = new Path(numVertices);
        boolean visited [] = new boolean [numVertices];
        int curVertex, closest = 0;
        double minDist;
        Arrays.fill(visited, false);

        shortest.path[0] = 0;
        visited[0] = true;

        // for each vertex
        for(int i = 1; i < numVertices; i++)
        {
            curVertex = shortest.path[i-1];
            minDist = Double.MAX_VALUE;
            // find the closest vertex
            for(int j = 0; j < numVertices; j++)
            {
                if(j != curVertex && !visited[j] && costMatrix[curVertex][j] < minDist)
                {
                    closest = j;
                    minDist = costMatrix[curVertex][j];
                }
            }

            // put the closest vertex as the next vertex in the path
            shortest.path[i] = closest;
            shortest.pathCost += minDist;
            visited[closest] = true;
        }
        shortest.pathCost += costMatrix[closest][0];
        return shortest;
    }

    // generate a random cost matrix
    static double [][] GenerateRandomCostMatrix(int numVertices, int maxCost)
    {
        double [][] costMatrix = new double [numVertices][numVertices];
        //for each vertex
        for(int v1 = 0; v1 < numVertices; v1++)
        {
            // for each other vertex
            for(int v2 = v1; v2 < numVertices; v2++)
            {
                // asign a random number up to maxCost to the edge between the two vertices
                costMatrix[v1][v2] = Math.random()*(double)maxCost;
                costMatrix[v2][v1] = costMatrix[v1][v2];
            }
        }
        return costMatrix;
    }

    // generate a random euclidean cost matrix
    static double [][] GenerateRandomEuclideanCostMatrix(int numVertices, int maxX, int maxY)
    {
        double costMatrix [][] = new double [numVertices][numVertices];
        double x [] = new double [numVertices];
        double y [] = new double [numVertices];

        //for each vertex generate a random x,y coordinate
        for(int v = 0; v < numVertices; v++)
        {
            x[v] = Math.random() * (double)maxX;
            y[v] = Math.random() * (double)maxY;
        }

        // for each vertex
        for(int v1 = 0; v1 < numVertices; v1++)
        {
            // for each other vertex
            for(int v2 = v1; v2 < numVertices; v2++)
            {
                // calculate the distance between the vertexes for the edge cost and store in cost matrix
                costMatrix[v1][v2] = distance(x[v1], y[v1], x[v2], y[v2]);
                costMatrix[v1][v2] = Double.parseDouble((new DecimalFormat("#.############")).format(costMatrix[v1][v2] ));
                costMatrix[v2][v1] = costMatrix[v1][v2];
            }
        }
        return costMatrix;
    }

    // generate a random circular graph cost matrix
    static double [][] GenerateRandomCircularGraphCostMatrix(int numVertices,int radius)
    {
        double costMatrix [][] = new double [numVertices][numVertices];
        double x [] = new double [numVertices];
        double y [] = new double [numVertices];
        boolean vertexAssigned [] = new boolean[numVertices];
        Arrays.fill(vertexAssigned, false);
        int randomVertex;
        double stepAngle = (2.0*Math.PI)/(double)numVertices;
        double stepDistance = 0;

        // geneate the initial x and y coordinate for the circular graph
        x[0] = (double)radius * Math.sin(0);
        y[0] = (double)radius * Math.cos(0);
        vertexAssigned[0] = true;
    //    System.out.println("Vertices and coordinates:");
    //    System.out.println("=========================");
    //    System.out.println("Vertex       x             y");
    //    System.out.printf(" 0     %11.6f    %11.6f\n", x[0], y[0]);
    //    System.out.print("Expected Path      :  0   ");

        // for each step
        for(int s = 1; s < numVertices; s++)
        {
            // get a random vertex that hasn't been visited
            do {
                randomVertex = (int)(Math.random() * numVertices);
            }while(vertexAssigned[randomVertex]);

            // assign the next coordinate on the circular graph to the random vertex
            vertexAssigned[randomVertex] = true;
            x[randomVertex] = (double)radius * Math.sin((double)s * stepAngle);
            y[randomVertex] = (double)radius * Math.cos((double)s * stepAngle);

    //      System.out.printf("%2d     %11.6f    %11.6f\n", randomVertex, x[randomVertex], y[randomVertex]);
            System.out.printf("%2d   ", randomVertex);
    //        if(s == 1)
    //            stepDistance = distance(x[0],y[0],x[randomVertex],y[randomVertex]);
        }
    //    System.out.printf("\nExpected Cost = %11.6f\n", stepDistance*numVertices);

        // for each vertex
        for(int v1 = 0; v1 < numVertices; v1++)
        {
            // for each other vertex
            for(int v2 = v1; v2 < numVertices; v2++)
            {
                // calculate the distance between the vertexes for the edge cost and store in cost matrix
                costMatrix[v1][v2] = distance(x[v1], y[v1], x[v2], y[v2]);
                costMatrix[v1][v2] = Double.parseDouble((new DecimalFormat("#.##########")).format(costMatrix[v1][v2] ));
                costMatrix[v2][v1] = costMatrix[v1][v2];
            }
        }
        return costMatrix;
    }

    // ant colonny algorithm for TSP
    static Path AntColonyTSP( double [][] costMatrix, int numAnts, double phermoneFactor, int maxUnchangedSteps)
    {
        int numVertices = costMatrix[0].length;
        Path bestPath = new Path(numVertices);
        double [][] phero = new double [numVertices][numVertices];
        double [][] newPhero = new double [numVertices][numVertices];
        bestPath.pathCost = Double.MAX_VALUE;
        int k,  h=0;
        double cumProb, edgeSelectionProbabilty;
        double pathCost, Q, totalA;
        int path [] = new int [numVertices];
        boolean visited [] = new boolean [numVertices];
        //for each time  step
        boolean changed;
        int unchangedCount = 0;
        int timestep = 0;
        double decayFactor = 0.5;

        // while the unchanged count is less than the maximum unchanged steps
        while(unchangedCount < maxUnchangedSteps)
        {
            changed = false;
            // for each ant
            for(int ant = 0; ant < numAnts; ant++)
            {
                Arrays.fill(visited, false);
                pathCost = 0;
                path[0] = 0;
                visited[0] = true;

                // for each step
                for(int step = 1; step < numVertices; step++)
                {
                    k = path[step - 1];
                    totalA = 0;
                    // for each vertex that hasn't been visited
                    for (h = 0; h < numVertices; h++) {
                        // add the attraction of that node from the current node to the total attraction
                        if (!visited[h] && h != k)
                            totalA += attraction(phero, costMatrix, k, h);
                    }

                    Q = Math.random();
                    cumProb = 0;
                    // for each vertex that hasn't been visited
                    for (h = 0; h < numVertices; h++) {
                        if (!visited[h] && h != k)
                        {
                            // calculate the edge probability for the current vertex
                            edgeSelectionProbabilty = attraction(phero, costMatrix, k, h) / totalA;
                            cumProb = cumProb + edgeSelectionProbabilty;
                            // if the random number is less than the cumulative probability then break and the current h, is the next vertex visited
                            if (Q < cumProb) {
                                break;
                            }
                        }
                    }
                    path[step] = h;
                    visited[h] = true;
                    pathCost += costMatrix[k][h];
                }
                pathCost = pathCost + costMatrix[h][0];

                // if the path taken by the ant is better than the best found so far, update the path and cost
                if(pathCost < bestPath.pathCost || (timestep == 0&& ant == 0))
                {
                    bestPath.path = path.clone();
                    bestPath.pathCost = pathCost;
                    changed = true;
                }
                // add the phermones that the ant laid
                for(int step = 0; step < numVertices; step++)
                {
                    k = path[step];
                    h = path[(step+1)%numVertices];
                    newPhero[k][h] += phermoneFactor/pathCost;
                }
            }
            // update the phermones from the ansts in this time step
            for(k = 0; k < numVertices; k++)
            {
                for(h=0; h< numVertices; h++)
                {
                    phero[k][h] *= decayFactor;
                    phero[k][h] += newPhero[k][h];
                }
            }

            //update the count depending on if the best path was updated
            if (!changed) {
                unchangedCount++;
            }
            else{
                unchangedCount=0;
            }
            timestep++;
        }
        return bestPath;
    }

    // calculate the attraction an ant has a certain edge
    static double attraction(double [][] phero, double [][] costMatrix, int k, int h)
    {
        return (1+phero[k][h])/costMatrix[k][h];
    }

    // calculate the distance between two coordinates
    static double distance (double x1, double y1, double x2, double y2)
    {
        return Math.pow((Math.pow((x1-x2),2) + Math.pow((y1-y2),2)), 0.5);
    }

    // run the full experiment
    static void runFullExperiment(String resultsFileName)
    {
        try {
            resultsFile = new FileWriter(ResultsFolderPath + resultsFileName);
            resultsWriter = new PrintWriter(resultsFile);
        } catch(Exception e) {
            System.out.println("*****!!!!!  Had a problem opening the results file "+ResultsFolderPath+resultsFileName);
            return;
        }

        // create stopwatch for timing individual trials
        ThreadCpuStopWatch TrialStopwatch = new ThreadCpuStopWatch();

        // print the labels in the results file
        resultsWriter.println("#InputSize       AverageTime  "  );
        resultsWriter.flush();

        //print a message to indicate the function and run that is currently running
        System.out.println(resultsFileName);

        //for each size of input we want to test: starting at MININPUTSIZE and doubling each iteration until reaching MAXINPUTSIZE
        for(int inputSize=MININPUTSIZE;inputSize<=MAXINPUTSIZE; inputSize+= 25 )
        {
            Path brutePath = new Path(inputSize);
            // Path greedyPath = new Path(inputSize);
            // Path antColPath = new Path(inputSize);

            //print progress message
            System.out.println("Running test for input size "+inputSize+" ... ");
            System.out.print("    Running trial batch...");

            // reset elapsed time for the batch to 0
            long batchElapsedTime = 0;

            // force garbage collection before each batch of trials run
            System.gc();

            // repeat for desired number of trials (for a specific size of input)...
            for (long trial = 0; trial < numberOfTrials; trial++)
            {
                double [][] costMatrix = GenerateRandomCostMatrix(inputSize, 150);

                // being timing
                TrialStopwatch.start();

                bruteForceTSP(costMatrix);
                // GreedyTSP(costMatrix);
                // AntColonyTSP(costMatrix, (int)((double)inputSize * 0.75),1, 75);


                //stop the timer and add to the total time elapsed for the batch of trials
                batchElapsedTime = batchElapsedTime + TrialStopwatch.elapsedTime();
            }
            // calculate the average time per trial in this batch
            double averageTimePerTrialInBatch = (double) batchElapsedTime / (double)numberOfTrials;

            // print data for this size of input
            resultsWriter.printf("%3d    %16.2f    \n",inputSize, averageTimePerTrialInBatch);
            resultsWriter.flush();
            System.out.println(" ....done.");
        }
    }
}
