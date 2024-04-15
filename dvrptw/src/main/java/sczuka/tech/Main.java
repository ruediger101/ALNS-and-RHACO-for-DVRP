package sczuka.tech;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.nio.file.Paths;

public class Main {
    public static void main(String[] args) throws FileNotFoundException {
        DataModel.testDataPathAndFile = "/home/ruediger/git/ALNS-and-RHACO-for-DVRP/cleanedupTestdata.csv";
        DataModel.resultOutputPath = "/home/ruediger/git/ALNS-and-RHACO-for-DVRP";

        File csvOutputFile = Paths.get(DataModel.resultOutputPath, "urgency.csv").toFile();
        try (PrintWriter pw = new PrintWriter(csvOutputFile)) {
            pw.println("UrgencyMean;UrgencyMedian;UrgencyStd");
            pw.flush();
            DataModel.loadData();
            for (int j = 10; j <= DataModel.noCustomers; j += 50) {
                DataModel.noImminentRequest = j;
                for (int k = 0; k < 6; k += 2) {
                    DataModel.urgencyIncrease = k / 10.0;
                    for (int i = 0; i < 5; i++) {
                        try {
                            DataModel.seed = i;
                            DataModel.loadData();
                            int[] urgency = DataModel.getUrgencyOfRequests();
                            double mean = DataModel.calculateMean(urgency);
                            double median = DataModel.calculateMedian(urgency);
                            double std = DataModel.calculateStandardDeviation(urgency);
                            pw.println(mean + ";" + median + ";" + std);
                            pw.flush();
                        } catch (Exception ex) {
                            ex.printStackTrace();
                        }
                    }
                }
            }
        }

        try {
            AntColonyOptimization.main(args);
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        try {
            AdaptiveLargeNeighborhoodSearch.main(args);
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

}
