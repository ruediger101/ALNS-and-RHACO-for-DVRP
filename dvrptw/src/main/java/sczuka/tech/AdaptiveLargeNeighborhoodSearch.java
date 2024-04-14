package sczuka.tech;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.nio.file.Paths;
import java.time.Duration;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Random;
import java.util.stream.Stream;

public class AdaptiveLargeNeighborhoodSearch {
    protected static DVRP currentSolution;
    protected static DVRP bestSolution;
    private final List<DestroyMethod> destroyMethods;
    private final List<RepairMethod> repairMethods;
    private final Random random;
    private double[] combinationWeights;
    private final double decayRate;
    private final int nonImprovingIerations;
    private final int innerIterations;
    private final int twoOptCount;
    private final double temperatureCoefficient = -0.05 / Math.log(0.5);
    private static int noNodesToDestroy = 40;

    public static void main(String[] args) throws FileNotFoundException {
        File csvOutputFile = Paths.get(DataModel.resultOutputPath, "ALNS.csv").toFile();
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
                            new AdaptiveLargeNeighborhoodSearch(500, 200, 200, 0.9);
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

    public AdaptiveLargeNeighborhoodSearch(int nonImprovingIerations, int innerIterations, int twoOptCount, double decayRate) {
        this.nonImprovingIerations = nonImprovingIerations;
        this.innerIterations = innerIterations;
        this.twoOptCount = twoOptCount;
        this.decayRate = decayRate;
        this.destroyMethods = new ArrayList<>();
        this.repairMethods = new ArrayList<>();
        this.random = new Random();
        addMethods();
        currentSolution = new DVRP(0);
        RepairMethod gi = new GreedyInsert();
        currentSolution = gi.apply(currentSolution);

        System.err.print("Epoche: 0");

        solve();

        RepairMethod kr3 = new KRegretInsert(3);
        for (int epochTime : DataModel.epochs) {
            System.err.print(", " + epochTime);
            currentSolution.increaseEpoche(epochTime);
            currentSolution = gi.apply(currentSolution);
            currentSolution = kr3.apply(currentSolution);
            solve();
        }

        System.err.println("\n\nBest Solution");
        System.err.println("\tNo Nodes: " + bestSolution.getTotalNodeNumber());
        System.err.println("\tFunktionswert: " + bestSolution.getFunctionValue());
        System.err.println("\tFahrzeuge in Verwendung: " + bestSolution.getVehicles().size());
        bestSolution.getVehicles()
                .forEach(v -> System.err.println("\tRoute: " + Stream.concat(v.getCompletedRoute().stream(), v.getRoute().stream()).toList()));
    }

    // Method to add destroy and repair methods and initialize weights afterward
    public void addMethods() {
        this.destroyMethods.add(new RandomRemoval());
        this.destroyMethods.add(new TabuRemoval());
        this.destroyMethods.add(new WorstCustomerRemoval());
        this.destroyMethods.add(new SimilarCustomerRemoval());
        this.destroyMethods.add(new TimeBasedRemoval());
        this.destroyMethods.add(new TimeWindowRemoval());

        this.repairMethods.add(new GreedyInsert());
        this.repairMethods.add(new RandomInsert());
        this.repairMethods.add(new KRegretInsert(2));
        this.repairMethods.add(new KRegretInsert(3));
        this.repairMethods.add(new GreedyInsert(0.25));
        this.repairMethods.add(new KRegretInsert(2, 0.25));
        this.repairMethods.add(new KRegretInsert(3, 0.25));
    }

    private void initializeWeights() {
        int totalCombinations = destroyMethods.size() * repairMethods.size();
        combinationWeights = new double[totalCombinations];
        Arrays.fill(combinationWeights, 1.0 / totalCombinations);
    }

    public void solve() {
        initializeWeights();
        if (combinationWeights == null || combinationWeights.length == 0) {
            throw new IllegalStateException("Methods and weights could not be initialized before solving.");
        }
        bestSolution = currentSolution;
        Map<Integer, double[]> scoreList = new HashMap<>();

        long iter = 0;
        for (int i = 1; i <= nonImprovingIerations; i++) {
            iter++;

            int combinationIndex = selectMethod(combinationWeights);
            int destroyIndex = combinationIndex / repairMethods.size();
            int repairIndex = combinationIndex % repairMethods.size();

            DVRP destroyed = destroyMethods.get(destroyIndex).apply(currentSolution);
            DVRP repaired = repairMethods.get(repairIndex).apply(destroyed);

            double score;
            if (isBetter(repaired, currentSolution)) {
                currentSolution = repaired.copy();
                if (isBetter(currentSolution, bestSolution)) {

                    i = 1;
                    score = 1;
                    bestSolution = currentSolution.copy();
                } else {
                    score = 0.4;
                }
            } else {
                double temp = temperatureCoefficient * currentSolution.getFunctionValue();
                double acceptanceProbability = Math.exp((currentSolution.getFunctionValue() - repaired.getFunctionValue()) / temp);
                if (random.nextDouble() <= acceptanceProbability) {
                    currentSolution = repaired.copy();
                    score = 0.25;
                } else {
                    score = 0;
                }

            }

            scoreList.compute(combinationIndex, (k, v) -> {
                if (v == null)
                    v = new double[] { 0, 0 };

                v[0] += score;
                v[1]++;
                return v;
            });

            if (iter % innerIterations == 0) {
                updateWeights(scoreList);
                scoreList.clear();
            }

            if (iter % twoOptCount == 0) {
                currentSolution = twoOpt(currentSolution);
                if (isBetter(currentSolution, bestSolution))
                    bestSolution = currentSolution.copy();
            }

            if (i % (nonImprovingIerations / 3) == 0)
                currentSolution = bestSolution.copy();
        }
        currentSolution = bestSolution.copy();
    }

    private int selectMethod(double[] weights) {
        double totalWeight = 0.0;
        for (double weight : weights) {
            totalWeight += weight;
        }
        double randomValue = random.nextDouble() * totalWeight;
        for (int i = 0; i < weights.length; i++) {
            randomValue -= weights[i];
            if (randomValue <= 0.0) {
                return i;
            }
        }
        return weights.length - 1; // Should not reach here
    }

    private void updateWeights(Map<Integer, double[]> scoreList) {
        for (Iterator<Entry<Integer, double[]>> it = scoreList.entrySet().iterator(); it.hasNext();) {
            Entry<Integer, double[]> next = it.next();

            combinationWeights[next.getKey()] = (1.0 - decayRate) * combinationWeights[next.getKey()] + decayRate * next.getValue()[0] / next.getValue()[1];
        }

    }

    protected boolean isBetter(DVRP s1, DVRP s2) {
        return s1.getFunctionValue() < s2.getFunctionValue();
    }

    protected abstract static class DestroyMethod {
        public abstract DVRP apply(DVRP solution);
    }

    protected abstract static class RepairMethod {
        public abstract DVRP apply(DVRP solution);
    }

    protected DVRP twoOpt(DVRP dvrp) {
        List<Integer> activeNodes = dvrp.getNodesInRoutes();
        if (activeNodes.size() < 2)
            return dvrp;

        DVRP bestImprovement = dvrp.copy();

        boolean improvementFound = true;
        while (improvementFound) {
            improvementFound = false;
            for (int node : activeNodes) {
                DVRP copy = dvrp.copy();
                if (copy.newRouteForCustomer(node) && isBetter(copy, bestImprovement)) {
                    bestImprovement = copy;
                    improvementFound = true;
                }
            }

            for (int node1 : activeNodes) {
                for (int node2 : activeNodes) {
                    DVRP copy = dvrp.copy();
                    if (copy.swapCustomers(node1, node2) && isBetter(copy, bestImprovement)) {
                        bestImprovement = copy;
                        improvementFound = true;
                    }
                }
            }
        }

        return bestImprovement;
    }

    protected class RandomRemoval extends DestroyMethod {
        public DVRP apply(DVRP s) {
            DVRP copy = s.copy();
            while (copy.getUnroutedCustomers().size() < noNodesToDestroy && !copy.getNodesInRoutes().isEmpty()) {
                copy.removeCustomerFromRoute(copy.getNodesInRoutes().get(random.nextInt(copy.getNodesInRoutes().size())));
            }
            return copy;
        }
    }

    protected class TabuRemoval extends DestroyMethod {
        private static List<Integer> tabuList = new ArrayList<>();
        private static int maxTabuListSize = 10;

        public DVRP apply(DVRP s) {
            DVRP copy = s.copy();
            List<Integer> removableNodes = new ArrayList<>(copy.getNodesInRoutes());
            removableNodes.removeAll(tabuList);
            while (copy.getUnroutedCustomers().size() < noNodesToDestroy && !removableNodes.isEmpty()) {
                int node = removableNodes.remove(random.nextInt(removableNodes.size()));
                copy.removeCustomerFromRoute(node);
                tabuList.add(node);
                if (tabuList.size() > maxTabuListSize) {
                    tabuList.remove(0);
                }
            }
            return copy;
        }
    }

    protected class WorstCustomerRemoval extends DestroyMethod {
        public DVRP apply(DVRP s) {
            DVRP worstCustomer = s.copy();
            while (worstCustomer.getUnroutedCustomers().size() < noNodesToDestroy && !worstCustomer.getNodesInRoutes().isEmpty()) {
                Iterator<Integer> it = worstCustomer.getNodesInRoutes().iterator();
                int node = it.next();
                DVRP worstCustomerInnerLoop = worstCustomer.copy();
                worstCustomerInnerLoop.removeCustomerFromRoute(node);

                while (it.hasNext()) {
                    node = it.next();
                    DVRP tempS = worstCustomer.copy();
                    tempS.removeCustomerFromRoute(node);

                    if (tempS.getFunctionValue(true) < worstCustomerInnerLoop.getFunctionValue(true)) {
                        worstCustomerInnerLoop = tempS;
                    }
                }
                worstCustomer = worstCustomerInnerLoop;
            }
            return worstCustomer;
        }
    }

    protected class SimilarCustomerRemoval extends DestroyMethod {
        private static double alpha = 1.0;
        private static double beta = 1.0;
        private static double gamma = 1.0;

        public DVRP apply(DVRP s) {
            DVRP copy = s.copy();
            if (copy.getUnroutedCustomers().size() >= noNodesToDestroy || copy.getNodesInRoutes().isEmpty())
                return copy;

            List<Integer> removedNodes = new ArrayList<>();
            int firstNode = copy.getNodesInRoutes().get(random.nextInt(copy.getNodesInRoutes().size()));
            copy.removeCustomerFromRoute(firstNode);
            removedNodes.add(firstNode);

            while (copy.getUnroutedCustomers().size() < noNodesToDestroy && !copy.getNodesInRoutes().isEmpty()) {
                int i = removedNodes.get(random.nextInt(removedNodes.size()));
                int mostSimilar = 0;
                double closestDistance = Double.MAX_VALUE;
                for (int k : copy.getNodesInRoutes()) {
                    double distance = alpha * DataModel.timeMatrix[i][k] + beta * Math.abs(DataModel.timeWindows[i][1] - DataModel.timeWindows[k][1])
                            + gamma * Math.abs(DataModel.timeWindows[i][0] - DataModel.timeWindows[k][0]);
                    if (distance < closestDistance) {
                        mostSimilar = k;
                        closestDistance = distance;
                    }
                }
                copy.removeCustomerFromRoute(mostSimilar);
            }

            return copy;
        }
    }

    protected class TimeBasedRemoval extends DestroyMethod {
        public DVRP apply(DVRP s) {
            DVRP copy = s.copy();
            int minTime = Arrays.asList(DataModel.timeWindows).stream().mapToInt(t -> t[0]).min().orElse(0);
            int maxTime = Arrays.asList(DataModel.timeWindows).stream().mapToInt(t -> t[1]).max().orElse(0);
            int randomTime = random.nextInt(maxTime - minTime) + minTime;

            List<Integer> list = new ArrayList<>(copy.getNodesInRoutes().stream().sorted((i, j) -> {
                int iDistance;
                int jDistance;

                if (randomTime < DataModel.timeWindows[i][0]) {
                    iDistance = DataModel.timeWindows[i][0] - randomTime;
                } else if (randomTime > DataModel.timeWindows[i][1]) {
                    iDistance = randomTime - DataModel.timeWindows[i][1];
                } else {
                    iDistance = 0;
                }

                if (randomTime < DataModel.timeWindows[j][0]) {
                    jDistance = DataModel.timeWindows[j][0] - randomTime;
                } else if (randomTime > DataModel.timeWindows[j][1]) {
                    jDistance = randomTime - DataModel.timeWindows[j][1];
                } else {
                    jDistance = 0;
                }

                return Integer.compare(iDistance, jDistance);
            }).toList());

            while (copy.getUnroutedCustomers().size() < noNodesToDestroy && !list.isEmpty()) {
                copy.removeCustomerFromRoute(list.remove(0));
            }

            return copy;

        }
    }

    protected class TimeWindowRemoval extends DestroyMethod {
        public DVRP apply(DVRP s) {
            DVRP result = s.copy();
            List<Integer> list = new ArrayList<>(result.getNodesInRoutes().stream().sorted((i, j) -> {
                int iDistance = DataModel.timeWindows[i][1] - DataModel.timeWindows[i][0];
                int jDistance = DataModel.timeWindows[j][1] - DataModel.timeWindows[j][0];

                return Integer.compare(iDistance, jDistance);
            }).toList());

            while (result.getUnroutedCustomers().size() < noNodesToDestroy && !list.isEmpty()) {
                result.removeCustomerFromRoute(list.remove(0));
            }

            return result;

        }
    }

    protected class GreedyInsert extends RepairMethod {
        private double eta;

        public GreedyInsert() {
            this(0);
        }

        public GreedyInsert(double eta) {
            this.eta = eta;
        }

        public DVRP apply(DVRP s) {
            DVRP result = s.copy();
            for (Iterator<Integer> it = result.getUnroutedCustomers().iterator(); it.hasNext();) {
                int node = it.next();
                DVRP tempResult = result.copy();
                it.remove();

                Vehicle nv = tempResult.getNewVehicle();
                nv.addToRoute(node);

                long bestValue;
                if (eta > 0)
                    bestValue = nv.getCost() + Math.round(eta * (random.nextDouble() * 2 - 1) * DataModel.maxTravelTime);
                else
                    bestValue = nv.getCost();

                for (int i = 0; i < result.getVehicles().size(); i++) {
                    Vehicle v = result.getVehicle(i);
                    for (int j = 0; j <= v.getRoute().size(); j++) {
                        Vehicle tempV = v.copy(null);
                        tempV.addToRoute(j, node);
                        long value;
                        if (eta > 0)
                            value = tempV.getCost() - v.getCost() + Math.round(eta * (random.nextDouble() * 2 - 1) * DataModel.maxTravelTime);
                        else
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
    }

    protected class RandomInsert extends RepairMethod {
        public DVRP apply(DVRP s) {
            DVRP result = s.copy();
            for (Iterator<Integer> it = result.getUnroutedCustomers().iterator(); it.hasNext();) {
                int node = it.next();

                List<int[]> validPositions = new ArrayList<>();
                if (result.isUnusedVehicleInDepot())
                    validPositions.add(new int[] { -1, 0 });

                for (int i = 0; i < result.getVehicles().size(); i++) {
                    Vehicle v = result.getVehicle(i);
                    for (int j = 0; j <= v.getRoute().size(); j++) {
                        if (v.copy(null).addToRoute(j, node)) {
                            validPositions.add(new int[] { i, j });
                        }
                    }
                }

                if (!validPositions.isEmpty()) {
                    it.remove();

                    int[] position = validPositions.get(random.nextInt(validPositions.size()));
                    if (position[0] == -1) {
                        result.getNewVehicle().addToRoute(node);
                    } else {
                        result.getVehicle(position[0]).addToRoute(position[1], node);
                    }
                }
            }
            return result;
        }
    }

    protected class KRegretInsert extends RepairMethod {
        private int k;
        private double eta;

        public KRegretInsert(int k) {
            this(k, 0);

        }

        public KRegretInsert(int k, double eta) {
            this.k = k;
            this.eta = eta;
        }

        public DVRP apply(DVRP s) {
            DVRP result = s.copy();
            while (!result.getUnroutedCustomers().isEmpty()) {

                DVRP bestResult = null;
                long maxRegret = -1;

                for (Iterator<Integer> it = result.getUnroutedCustomers().iterator(); it.hasNext();) {
                    int node = it.next();

                    Queue<Long> costIncrease = new PriorityQueue<>();
                    DVRP tempResult = result.copy();
                    Vehicle dv = tempResult.getNewVehicle();
                    dv.addToRoute(node);
                    if (eta > 0)
                        costIncrease.add(dv.getCost() + Math.round(eta * (random.nextDouble() * 2 - 1) * DataModel.maxTravelTime));
                    else
                        costIncrease.add(dv.getCost());

                    for (int i = 0; i < result.getVehicles().size(); i++) {
                        Vehicle v = result.getVehicle(i);
                        for (int j = 0; j <= v.getRoute().size(); j++) {
                            Vehicle tempV = v.copy(null);
                            tempV.addToRoute(j, node);
                            long cost;
                            if (eta > 0)
                                cost = tempV.getCost() - v.getCost() + Math.round(eta * (random.nextDouble() * 2 - 1) * DataModel.maxTravelTime);
                            else
                                cost = tempV.getCost() - v.getCost();

                            if (cost < costIncrease.peek()) {
                                tempResult = result.copy();
                                tempResult.getVehicle(i).addToRoute(j, node);
                            }
                            costIncrease.add(cost);
                        }
                    }

                    long regret;
                    if (costIncrease.size() >= k) {
                        long lowestValue = costIncrease.poll();
                        for (int i = 2; i < k; i++) {
                            costIncrease.poll();
                        }
                        regret = costIncrease.poll() - lowestValue;
                    } else {
                        regret = Long.MAX_VALUE;
                    }

                    if (regret > maxRegret) {
                        maxRegret = regret;
                        bestResult = tempResult;
                    }

                }
                result = bestResult;
            }
            return result;
        }
    }
}
