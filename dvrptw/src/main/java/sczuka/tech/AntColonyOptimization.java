package sczuka.tech;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.nio.file.Paths;
import java.time.Duration;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class AntColonyOptimization {
    private double[][] pheromones;
    private static final double initalPheromoneLevel = 0.1;
    private double pheromoneVotality = 1;
    private static final double minVotality = 0.2;
    private static double Q = 100;

    private static final double alpha = 1;
    private static final double beta = 3;
    private static final double gamma = 2;
    private static final double q_0 = 0.5;
    private static DVRP bestSolution = null;

    private static final Random rand = new Random(System.nanoTime());

    public static void main(String[] args) throws FileNotFoundException {
        File csvOutputFile = Paths.get(DataModel.resultOutputPath, "RHACO.csv").toFile();
        try (PrintWriter pw = new PrintWriter(csvOutputFile)) {
            pw.println("Imminent Request;dod;edod;dynamic;Funktionswert;Fahrzeuge");
            pw.flush();
            DataModel.loadData();
            for (int j = 10; j <= DataModel.noCustomers; j += 50) {
                DataModel.noImminentRequest = j;

                for (int k = 0; k < 6; k += 2) {
                    DataModel.urgencyIncrease = k / 10.0;
                    for (int i = 0; i < 5; i++) {
                        try {
                            System.err.println("###################################################################");
                            DataModel.seed = i;
                            DataModel.loadData();

                            System.err.println("Imminent Requests: " + j);
                            System.err.println("Seed: " + i);
                            System.err.println("Urgency Increase Factor: " + k);

                            Instant start = Instant.now();
                            new AntColonyOptimization();
                            Instant end = Instant.now();
                            System.err.println("\n\nDauer: " + Duration.between(start, end).toSeconds() + "s");
                            pw.println(DataModel.noImminentRequest + ";" + DataModel.dod + ";" + DataModel.edod + ";" + DataModel.dynamic + ";"
                                    + bestSolution.getFunctionValue() + ";" + bestSolution.getVehicles().size());
                            pw.flush();
                        } catch (Exception ex) {
                            ex.printStackTrace();
                        }
                    }
                }
            }
        }
    }

    private void initalizePheromones() {
        int size = DataModel.timeMatrix.length;
        pheromones = new double[size][size];
        IntStream.range(0, size).forEach(i -> IntStream.range(0, size).forEach(j -> pheromones[i][j] = initalPheromoneLevel));
    }

    private void updatePheromones(DVRP solution) {
        double deltaPheromones = Q / solution.getFunctionValue();

        solution.getVehicles().forEach(v -> {
            List<Integer> route = v.getRoute();
            if (!route.isEmpty()) {
                pheromones[v.getStartPostion()][route.get(0)] = (1.0 - pheromoneVotality) * pheromones[v.getStartPostion()][route.get(0)]
                        + pheromoneVotality * deltaPheromones;
                pheromones[route.get(route.size() - 1)][0] = (1.0 - pheromoneVotality) * pheromones[route.get(route.size() - 1)][0]
                        + pheromoneVotality * deltaPheromones;

                for (int i = 0; i < route.size() - 1; i++) {
                    pheromones[route.get(i)][route.get(i + 1)] = (1.0 - pheromoneVotality) * pheromones[route.get(i)][route.get(i + 1)]
                            + pheromoneVotality * deltaPheromones;
                }
            }
        });

    }

    public AntColonyOptimization() {
        int maxIterations = 5000;
        DataModel.loadData();
        initalizePheromones();

        System.err.print("Epoche: 0");

        DVRP s = new DVRP(0);
        s = HybridAntColonyOptimization(s, maxIterations);

        for (int epochTime : DataModel.epochs) {
            System.err.print(", " + epochTime);
            s.increaseEpoche(epochTime);
            List<double[][]> diversePheromones = pheromonDiversification(s, 10);
            ensembleMethod(s, diversePheromones, 10, 5);
            s = HybridAntColonyOptimization(s, maxIterations);
        }

        if (bestSolution != null) {
            System.err.println("\n\nBest Solution");
            System.err.println("\tNo Nodes: " + bestSolution.getTotalNodeNumber());
            System.err.println("\tFunktionswert: " + bestSolution.getFunctionValue());
            System.err.println("\tFahrzeuge in Verwendung: " + bestSolution.getVehicles().size());
            bestSolution.getVehicles()
                    .forEach(v -> System.err.println("\tRoute: " + Stream.concat(v.getCompletedRoute().stream(), v.getRoute().stream()).toList()));
        } else {
            System.err.println("No Solution found!");
        }
    }

    private void ensembleMethod(DVRP s, List<double[][]> pheromoneMatrixList, int matingListSize, int matingListSizeOfBest) {
        List<DVRP> foundSolutions = new ArrayList<>();

        DVRP noNewNodes = s.copy();
        noNewNodes.getUnroutedCustomers().clear();

        List<Integer> allNodes = new ArrayList<>(s.getNodesInRoutes());
        allNodes.addAll(s.getUnroutedCustomers());

        for (double[][] matrix : pheromoneMatrixList) {
            for (int i = 0; i < 4; i++) {
                pheromones = copyMatrix(matrix);

                DVRP newSolution = HybridAntColonyOptimization(noNewNodes.copy(), 1);

                if (!s.getNodesInRoutes().containsAll(newSolution.getNodesInRoutes()) || !newSolution.getNodesInRoutes().containsAll(s.getNodesInRoutes()))
                    System.err.println("HybridAntColonyOptimization Shit");

                newSolution.getUnroutedCustomers().addAll(s.getUnroutedCustomers());

                List<Integer> allNodes2 = new ArrayList<>(newSolution.getNodesInRoutes());
                allNodes2.addAll(newSolution.getUnroutedCustomers());
                if (!allNodes.containsAll(allNodes2) || !allNodes2.containsAll(allNodes))
                    System.err.println("Add Unrouted Shit");

                newSolution = greedyInsert(newSolution);

                if (!allNodes.containsAll(newSolution.getNodesInRoutes()) || !newSolution.getNodesInRoutes().containsAll(allNodes))
                    System.err.println("Greedy Shit");

                newSolution = twoOpt(newSolution);
                if (newSolution.isValid()) {
                    foundSolutions.add(newSolution.copy());
                }
            }
        }

        for (int k = 0; k < 20; k++) {
            List<DVRP> matingSolutions = new ArrayList<>(foundSolutions);
            for (int i = 0; i < foundSolutions.size() - 1; i++) {
                for (int j = i; j < foundSolutions.size(); j++) {
                    matingSolutions.addAll(mating(foundSolutions.get(i), foundSolutions.get(j)));
                }
            }
            foundSolutions.addAll(matingSolutions);
            foundSolutions.sort((a, b) -> Double.compare(a.getFunctionValue(), b.getFunctionValue()));
            if (foundSolutions.size() > matingListSize) {
                List<DVRP> temp = foundSolutions.subList(30, foundSolutions.size());
                Collections.shuffle(temp);

                foundSolutions = foundSolutions.subList(0, 30);
                foundSolutions.addAll(temp.subList(0, matingListSize - matingListSizeOfBest));
            }
        }

        List<DVRP> result = new ArrayList<>(foundSolutions);

        for (DVRP dvrp : foundSolutions) {
            DVRP newSolution = mutateSolution(dvrp);
            if (newSolution.isValid())
                result.add(newSolution);
        }

        result.sort((a, b) -> Double.compare(a.getFunctionValue(), b.getFunctionValue()));

        result = result.subList(0, matingListSize);

        initalizePheromones();
        for (DVRP dvrp : result) {
            updatePheromones(dvrp);
        }
    }

    private List<DVRP> mating(DVRP dvrp1, DVRP dvrp2) {
        DVRP a = dvrp1.copy();
        List<Vehicle> vaList = a.getVehicles().stream().filter(v -> !v.getRoute().isEmpty()).toList();
        Vehicle va = vaList.get(rand.nextInt(vaList.size()));
        int va1 = rand.nextInt(va.getRoute().size());
        int va2 = rand.nextInt(va.getRoute().size());
        List<Integer> nodesA = new ArrayList<>();

        if (va1 == va2)
            nodesA.add(va.getRoute().get(va1));
        else if (va1 < va2)
            nodesA.addAll(va.getRoute().subList(va1, va2));
        else
            nodesA.addAll(va.getRoute().subList(va2, va1));

        DVRP b = dvrp2.copy();
        List<Vehicle> vbList = b.getVehicles().stream().filter(v -> !v.getRoute().isEmpty()).toList();
        Vehicle vb = vbList.get(rand.nextInt(vbList.size()));
        int vb1 = rand.nextInt(vb.getRoute().size());
        int vb2 = rand.nextInt(vb.getRoute().size());
        List<Integer> nodesB = new ArrayList<>();

        if (vb1 == vb2)
            nodesB.add(vb.getRoute().get(vb1));
        else if (vb1 < vb2)
            nodesB.addAll(vb.getRoute().subList(vb1, vb2));
        else
            nodesB.addAll(vb.getRoute().subList(vb2, vb1));

        List<DVRP> result = new ArrayList<>();

        a.removeCustomersFromRoute(nodesB);
        a = greedyInsert(a, nodesB);
        if (a.isValid())
            result.add(a);

        b.removeCustomersFromRoute(nodesA);
        b = greedyInsert(b, nodesA);
        if (b.isValid())
            result.add(b);

        return result;
    }

    private List<double[][]> pheromonDiversification(DVRP s, int n_s) {
        int unvisited = s.getUnroutedCustomers().size() + s.getNodesInRoutes().size();
        double d_s = ((double) s.getUnroutedCustomers().size()) / unvisited;
        int n_m = (int) Math.ceil(d_s * (unvisited * (unvisited - 1)) / 2.0);

        double[][][] distanceMatrix = new double[s.getNodesInRoutes().size()][s.getNodesInRoutes().size()][2];

        double maxDistance = 0;

        for (int i = 0; i < s.getNodesInRoutes().size() - 1; i++) {
            int nodeA = s.getNodesInRoutes().get(i);
            for (int j = i + 1; j < s.getNodesInRoutes().size() - 1; j++) {
                int nodeB = s.getNodesInRoutes().get(j);
                int nearestNode = -1;
                double smallestDistance = Double.MAX_VALUE;
                for (int nodeC : s.getUnroutedCustomers()) {
                    int c = DataModel.timeMatrix[nodeA][nodeB];
                    int b = DataModel.timeMatrix[nodeA][nodeC];
                    int a = DataModel.timeMatrix[nodeB][nodeC];

                    double[] angles = TriangleCalculations.calculateAngles(a, b, c);
                    double distance;
                    if (angles[0] >= 90.0)
                        distance = b;
                    else if (angles[1] >= 90.0)
                        distance = a;
                    else {
                        double area = TriangleCalculations.calculateArea(a, b, c);
                        distance = TriangleCalculations.calculateHeight(area, a);
                    }

                    if (distance < smallestDistance) {
                        smallestDistance = distance;
                        nearestNode = nodeC;
                    }
                }
                distanceMatrix[i][j][0] = smallestDistance;
                distanceMatrix[i][j][1] = nearestNode;
                distanceMatrix[j][i][0] = smallestDistance;
                distanceMatrix[j][i][1] = nearestNode;
                if (smallestDistance > maxDistance)
                    maxDistance = smallestDistance;
            }
        }

        List<double[][]> modifiedPheromoneMatices = new ArrayList<>();

        List<int[]> customerPairs = new ArrayList<>();
        for (int i = 0; i < s.getNodesInRoutes().size() - 1; i++) {
            for (int j = i + 1; j < s.getNodesInRoutes().size() - 1; j++) {
                if (distanceMatrix[i][j][0] > 0.001)
                    customerPairs.add(new int[] { i, j });
            }
        }

        double meanPheromones = meanOfMatrix(pheromones);
        Map<int[], Integer> shortestArcs = new HashMap<>();
        for (int m = 0; m < n_s; m++) {
            double[][] modifiedPheromones = copyMatrix(pheromones);
            List<int[]> localCustomerPairs = new ArrayList<>(customerPairs);
            Set<int[]> selectedPairs = new HashSet<>();

            if (customerPairs.size() > n_m) {
                while (selectedPairs.size() < n_m) {
                    Collections.shuffle(localCustomerPairs);

                    for (Iterator<int[]> it = localCustomerPairs.iterator(); it.hasNext();) {
                        int[] pair = it.next();
                        double p = rand.nextDouble();
                        if (p <= distanceMatrix[pair[0]][pair[1]][0] / maxDistance) {
                            selectedPairs.add(pair);
                            it.remove();
                        }
                    }
                }
            } else {
                selectedPairs.addAll(customerPairs);
            }

            for (int[] pair : selectedPairs) {
                if (modifiedPheromones[pair[0]][pair[1]] > meanPheromones) {
                    int nodeA = s.getNodesInRoutes().get(pair[0]);
                    int nodeB = s.getNodesInRoutes().get(pair[1]);

                    shortestArcs.computeIfAbsent(pair, (v) -> {
                        int newNode = (int) (distanceMatrix[pair[0]][pair[1]][1] + 0.5);
                        int abc = DataModel.timeMatrix[nodeA][nodeB] + DataModel.timeMatrix[nodeB][newNode];
                        int cab = DataModel.timeMatrix[newNode][nodeA] + DataModel.timeMatrix[nodeA][nodeB];
                        int acb = DataModel.timeMatrix[nodeA][newNode] + DataModel.timeMatrix[newNode][nodeB];
                        int shortestArc = Math.min(abc, Math.min(cab, acb));
                        return shortestArc;
                    });

                    double mt = shortestArcs.get(pair) / (DataModel.timeMatrix[nodeA][nodeB] + Double.MIN_VALUE);
                    modifiedPheromones[nodeA][nodeB] *= mt;
                    modifiedPheromones[nodeB][nodeA] *= mt;
                }
            }
            modifiedPheromoneMatices.add(modifiedPheromones);
        }
        return modifiedPheromoneMatices;
    }

    private double[][] copyMatrix(double[][] p) {
        double[][] copy = new double[p.length][p[0].length];

        for (int i = 0; i < p.length; i++) {
            for (int j = 0; j < p[0].length; j++) {
                copy[i][j] = p[i][j];
            }
        }
        return copy;
    }

    private double meanOfMatrix(double[][] p) {
        double sum = 0;
        int count = 0;
        for (int i = 0; i < p.length - 1; i++) {
            for (int j = i + 1; j < p[0].length; j++) {
                sum += p[i][j];
                count++;
            }
        }
        return sum / count;
    }

    public DVRP HybridAntColonyOptimization(DVRP s, int maxIterations) {
        double bestFunctionValue = Double.MAX_VALUE;
        bestSolution = null;
        DVRP emptyRoutes = s.copy();
        emptyRoutes.removeUnvisited();

        for (int i = 0; i < maxIterations; i++) {
            DVRP tempS = emptyRoutes.copy();

            Collections.shuffle(tempS.getVehicles());

            for (Vehicle v : tempS.getVehicles()) {
                while (!tempS.getUnroutedCustomers().isEmpty()) {
                    int j = getNextNode(tempS.getUnroutedCustomers(), v.getPosition());
                    if (v.isValidRoute(j)) {
                        v.addToRoute(j);
                    } else {
                        break;
                    }
                }
            }
            while (!tempS.getUnroutedCustomers().isEmpty() && tempS.isUnusedVehicleInDepot()) {
                Vehicle v = tempS.getNewVehicle();
                while (!tempS.getUnroutedCustomers().isEmpty()) {
                    int j = getNextNode(tempS.getUnroutedCustomers(), v.getPosition());
                    if (v.isValidRoute(j)) {
                        v.addToRoute(j);
                    } else {
                        break;
                    }
                }
                if (v.getRoute().isEmpty())
                    tempS.getVehicles().remove(v);
            }

            if (!tempS.getUnroutedCustomers().isEmpty()) {
                continue;
            }

            if (tempS.getFunctionValue() < bestFunctionValue) {
                bestSolution = tempS.copy();
                bestFunctionValue = tempS.getFunctionValue();
            }

            tempS = mutateSolution(tempS);
            if (tempS != null && tempS.isValid()) {
                if (tempS.getFunctionValue() < bestFunctionValue) {
                    bestSolution = tempS.copy();
                    bestFunctionValue = tempS.getFunctionValue();
                }
                updatePheromones(tempS);
                pheromoneVotality = Math.max(0.9 * pheromoneVotality, minVotality);
            }

        }
        return bestSolution;
    }

    private DVRP mutateSolution(DVRP solution) {
        DVRP copy = solution.copy();

        // Swap two customers
        // Include depot in random numbers to allow creation of new routes

        List<Integer> activeNodes = copy.getNodesInRoutes();
        if (activeNodes.size() < 2)
            return null;

        activeNodes.add(0);
        int first = activeNodes.remove(rand.nextInt(activeNodes.size()));
        int second = activeNodes.remove(rand.nextInt(activeNodes.size()));

        // Create new route if one of the selection is the depot
        if (first == 0 || second == 0) {
            copy.newRouteForCustomer(Math.max(first, second));
        } else {
            copy.swapCustomers(first, second);
        }

        // Insert Customer at another location
        activeNodes = copy.getNodesInRoutes();
        int nodeToInsert = activeNodes.remove(rand.nextInt(activeNodes.size()));
        int insertBeforeNode = activeNodes.remove(rand.nextInt(activeNodes.size()));

        copy.addCustomerBeforeOther(nodeToInsert, insertBeforeNode);

        return copy;
    }

    private int getNextNode(Set<Integer> unvisitedCustomer, int i) {
        double q = rand.nextDouble();
        Map<Integer, Double> args = new HashMap<>();
        for (int j : unvisitedCustomer) {
            double arg = Math.pow(pheromones[i][j], alpha) * Math.pow(1.0 / DataModel.timeMatrix[i][j], beta)
                    * Math.pow(1.0 / (DataModel.timeWindows[j][1] - DataModel.timeWindows[j][0]), gamma);
            args.put(j, arg);
        }

        if (q <= q_0) {
            return args.entrySet().stream().sorted((a, b) -> Double.compare(b.getValue(), a.getValue())).map(Entry::getKey).findFirst().orElse(-1);
        } else {
            double argsSum = args.values().stream().reduce(0.0, Double::sum);
            double r = rand.nextDouble();
            double summedProbabilities = 0;
            for (Entry<Integer, Double> e : args.entrySet()) {
                summedProbabilities += (e.getValue() / argsSum);
                if (r <= summedProbabilities)
                    return e.getKey();
            }

            return args.entrySet().stream().toList().get(args.size() - 1).getKey();
        }
    }

    public DVRP greedyInsert(DVRP s) {
        DVRP result = s.copy();
        for (Iterator<Integer> it = result.getUnroutedCustomers().iterator(); it.hasNext();) {
            int node = it.next();
            DVRP tempResult = result.copy();
            it.remove();

            Vehicle nv = tempResult.getNewVehicle();
            nv.addToRoute(node);

            long bestValue;
            bestValue = nv.getCost();

            for (int i = 0; i < result.getVehicles().size(); i++) {
                Vehicle v = result.getVehicle(i);
                for (int j = 0; j <= v.getRoute().size(); j++) {
                    Vehicle tempV = v.copy(null);
                    tempV.addToRoute(j, node);
                    long value;
                    value = tempV.getCost() - v.getCost();

                    if (value < bestValue) {
                        bestValue = value;
                        tempResult = result.copy();
                        tempResult.getVehicle(i).addToRoute(j, node);
                    }
                }
            }

            result = tempResult;

        }
        return result;
    }

    public DVRP greedyInsert(DVRP s, List<Integer> customers) {
        DVRP result = s.copy();
        Vehicle nv = result.getNewVehicle();
        nv.addToRoute(customers);

        long bestValue;
        bestValue = nv.getCost();

        for (int i = 0; i < s.getVehicles().size(); i++) {
            Vehicle v = s.getVehicle(i);
            for (int j = 0; j <= v.getRoute().size(); j++) {
                Vehicle tempV = v.copy(null);
                tempV.addToRoute(j, customers);
                long value;
                value = tempV.getCost() - v.getCost();

                if (value < bestValue) {
                    bestValue = value;
                    result = s.copy();
                    result.getVehicle(i).addToRoute(j, customers);
                }
            }
        }
        return result;
    }

    protected DVRP twoOpt(DVRP dvrp) {
        List<Integer> activeNodes = dvrp.getNodesInRoutes();
        if (activeNodes.size() < 2)
            return dvrp;

        DVRP bestImprovement = dvrp.copy();

        boolean improvementFound = true;
        while (improvementFound) {
            improvementFound = false;
            for (Integer node : activeNodes) {
                DVRP copy = dvrp.copy();
                if (copy.newRouteForCustomer(node) && copy.getFunctionValue() < bestImprovement.getFunctionValue()) {
                    bestImprovement = copy;
                    improvementFound = true;
                }
            }

            for (int node1 : activeNodes) {
                for (int node2 : activeNodes) {
                    DVRP copy = dvrp.copy();
                    if (copy.swapCustomers(node1, node2) && copy.getFunctionValue() < bestImprovement.getFunctionValue()) {
                        bestImprovement = copy;
                        improvementFound = true;
                    }
                }
            }
        }

        return bestImprovement;
    }

}
