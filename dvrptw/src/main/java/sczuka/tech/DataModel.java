package sczuka.tech;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.NumberFormat;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Random;
import java.util.stream.IntStream;

class DataModel {
    public static void loadData() {

        try {
            List<String> allLines = Files.readAllLines(Paths.get(testDataPathAndFile));

            allLines.removeIf(s -> s.isBlank());

            demands = new int[allLines.size()];
            serviceTimes = new int[allLines.size()];
            timeWindows = new int[allLines.size()][2];
            timeMatrix = new int[allLines.size()][allLines.size()];
            requestArrivalTime = new int[allLines.size()];

            NumberFormat format = NumberFormat.getInstance(Locale.GERMANY);

            for (int i = 0; i < allLines.size(); i++) {
                String[] columns = allLines.get(i).split(";");
                serviceTimes[i] = format.parse(columns[3]).intValue();
                demands[i] = format.parse(columns[4]).intValue();
                timeWindows[i][0] = format.parse(columns[5]).intValue();
                timeWindows[i][1] = format.parse(columns[6]).intValue();

                for (int j = 0; j < allLines.size(); j++) {
                    String[] columns1 = allLines.get(i).split(";");
                    String[] columns2 = allLines.get(j).split(";");
                    double x1 = format.parse(columns1[1]).doubleValue();
                    double y1 = format.parse(columns1[2]).doubleValue();
                    double x2 = format.parse(columns2[1]).doubleValue();
                    double y2 = format.parse(columns2[2]).doubleValue();

                    int travelTime = (int) Math.round(Math.sqrt(Math.pow(Math.abs(x1 - x2), 2) + Math.pow(Math.abs(y1 - y2), 2)));
                    maxTravelTime = Math.max(maxTravelTime, travelTime);
                    timeMatrix[i][j] = travelTime;

                }
            }

            functionPenaltyValue = Arrays.stream(DataModel.timeMatrix).mapToLong(i -> Arrays.stream(i).mapToLong(j -> j).max().orElse(0)).sum()
                    + vehicleNumber * vehicleFixCost;

            requestArrivalTime[0] = 0; // Depot always zero, only in list to have identical indices across lists

            List<Integer> customerIDs = new ArrayList<>(IntStream.range(1, allLines.size()).boxed().toList());
            Collections.shuffle(customerIDs);

            Iterator<Integer> it = customerIDs.iterator();
            for (int i = 0; i < noImminentRequest && it.hasNext(); i++) {
                requestArrivalTime[it.next()] = 0;
                it.remove();
            }

            List<Integer> dynamicCustomers = new ArrayList<>(customerIDs);

            Random r = new Random(seed);
            while (it.hasNext()) {
                int i = it.next();
                it.remove();
                int possibleRequestWindow = timeWindows[i][1] - timeMatrix[0][i];
                requestArrivalTime[i] = r.nextInt((int) Math.round(urgencyIncrease * possibleRequestWindow), possibleRequestWindow) + 1;
            }

            double edodSumValue = 0;
            endOfPlanning = Arrays.asList(timeWindows).stream().mapToInt(t -> t[1]).max().orElse(0);

            epochs = new ArrayList<>();
            List<Integer> requestList = new ArrayList<>(Arrays.stream(requestArrivalTime).boxed().toList());
            requestList.removeIf(req -> req == 0);

            int noNewVehicles = 0;
            for (int t = 1; t <= endOfPlanning; t++) {
                for (it = requestList.iterator(); it.hasNext();) {
                    if (it.next() == t) {
                        it.remove();
                        noNewVehicles++;
                    }
                }
                if (noNewVehicles >= vehiclesInEpoch) {
                    epochs.add(t);
                    noNewVehicles = 0;
                }
            }

            for (it = dynamicCustomers.iterator(); it.hasNext();) {
                int node = it.next();
                int urgency = timeWindows[node][1] - requestArrivalTime[node];
                edodSumValue += (1.0 - ((double) urgency) / endOfPlanning);
            }

            List<Integer> sortedArrivalTime = dynamicCustomers.stream().map(cust -> requestArrivalTime[cust]).sorted().toList();
            List<Integer> deltaArrivalTimes = new ArrayList<>();
            for (int i = 0; i < sortedArrivalTime.size() - 1; i++) {
                deltaArrivalTimes.add(sortedArrivalTime.get(i + 1) - sortedArrivalTime.get(i));
            }

            Double perfectInterval = ((double) endOfPlanning) / deltaArrivalTimes.size();
            List<Double> sigma = new ArrayList<>();
            if (deltaArrivalTimes.get(0) < perfectInterval)
                sigma.add(perfectInterval - deltaArrivalTimes.get(0));
            else
                sigma.add(0.0);

            List<Double> maxSigma = new ArrayList<>(sigma);

            for (int i = 1; i < deltaArrivalTimes.size(); i++) {
                if (deltaArrivalTimes.get(i) < perfectInterval) {
                    double diff = perfectInterval - deltaArrivalTimes.get(i);
                    double penalty = diff * sigma.get(i - 1) / perfectInterval;
                    sigma.add(diff + penalty);
                    maxSigma.add(perfectInterval + penalty);
                } else {
                    sigma.add(0.0);
                    maxSigma.add(0.0);
                }
            }

            noCustomers = allLines.size() - 1;
            dod = ((double) noCustomers - noImminentRequest) / noCustomers;
            edod = edodSumValue / noCustomers;
            dynamic = 1.0 - sigma.stream().mapToDouble(Double::doubleValue).sum() / maxSigma.stream().mapToDouble(Double::doubleValue).sum();

        } catch (IOException | SecurityException | ParseException e) {
            e.printStackTrace();
        }
    }

    public static double calculateStandardDeviation(int[] array) {
        double mean = calculateMean(array);

        // calculate the standard deviation
        double standardDeviation = 0.0;
        for (double num : array) {
            standardDeviation += Math.pow(num - mean, 2);
        }

        return Math.sqrt(standardDeviation / array.length);
    }

    public static double calculateMean(int[] array) {
        // get the sum of array
        double sum = 0.0;
        for (double i : array) {
            sum += i;
        }

        // get the mean of array
        int length = array.length;
        return sum / length;
    }

    public static double calculateMedian(int[] array) {
        int length = array.length;
        int[] copy = Arrays.copyOf(array, length);

        Arrays.sort(copy);
        // get the sum of array
        int floorOfHalfLength = length / 2;
        if (length % 2 == 0) {
            // get the mean of the two value in the middle of the sorted array
            return (copy[floorOfHalfLength] + copy[floorOfHalfLength - 1]) / 2.0;

        } else {
            // get the value in the middle of the sorted array
            return copy[floorOfHalfLength];
        }
    }

    public static int[] getUrgencyOfRequests() {
        int[] urgency = new int[DataModel.requestArrivalTime.length - 1];
        for (int x = 1; x < DataModel.requestArrivalTime.length; x++) {
            urgency[x - 1] = DataModel.timeWindows[x][1] - DataModel.requestArrivalTime[x];
        }
        return urgency;
    }

    public static int[] requestArrivalTime;
    public static int[][] timeMatrix;
    public static int[][] timeWindows;
    // [START demands_capacities]
    public static int[] demands;
    public static int[] serviceTimes;
    public static int vehicleCapacities = 250;
    public static int vehicleFixCost = 5;

    public static int vehicleNumber = 40;
    public static int depot = 0;
    public static long maxTravelTime = 0;
    public static boolean allowInvalidRoutes = true;
    public static double lambda = 0.1;
    public static long functionPenaltyValue;
    public static long seed = 0;
    public static int noImminentRequest = 0;
    public static int vehiclesInEpoch = 5;
    public static double urgencyIncrease = 0;
    public static int endOfPlanning = 0;
    public static int noCustomers = 0;
    public static double dod = 0;
    public static double edod = 0;
    public static double dynamic = 0;
    public static List<Integer> epochs = new ArrayList<>();
    public static String testDataPathAndFile;
    public static String resultOutputPath;

}
