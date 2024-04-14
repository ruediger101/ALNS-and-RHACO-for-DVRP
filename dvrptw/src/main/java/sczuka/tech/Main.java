package sczuka.tech;

import java.io.FileNotFoundException;

public class Main {
    public static void main(String[] args) throws FileNotFoundException {
        DataModel.testDataPathAndFile = "/home/ruediger/git/ALNS-and-RHACO-for-DVRP/cleanupTestdata.csv";
        DataModel.resultOutputPath = "/home/ruediger/git/ALNS-and-RHACO-for-DVRP";
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
