package sczuka.tech;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

class DVRP {
    private List<Vehicle> vehicles = new ArrayList<>();
    private Set<Integer> unroutedCustomers = new HashSet<>();
    private int noVehicles = DataModel.vehicleNumber;
    private int vehicleCapacity = DataModel.vehicleCapacities;
    private int vehicleFixCost = DataModel.vehicleFixCost;
    private int epochStartTime = 0;

    DVRP() {
    }

    DVRP(int startNewEpoch) {
        addCustomerRequestInInterval(-1, startNewEpoch);
    }

    private void addCustomerRequestInInterval(int startLastEpoch, int startNewEpoch) {
        for (int i = 1; i < DataModel.requestArrivalTime.length; i++) {
            if (DataModel.requestArrivalTime[i] > startLastEpoch && DataModel.requestArrivalTime[i] <= startNewEpoch)
                unroutedCustomers.add(i);
        }
        epochStartTime = startNewEpoch;
    }

    public DVRP copy() {
        DVRP copy = new DVRP();
        copy.vehicles = new ArrayList<>();
        for (Vehicle v : vehicles) {
            copy.vehicles.add(v.copy(copy));
        }
        copy.unroutedCustomers = new HashSet<>(unroutedCustomers);
        copy.noVehicles = noVehicles;
        copy.vehicleCapacity = vehicleCapacity;
        copy.vehicleFixCost = vehicleFixCost;
        copy.epochStartTime = epochStartTime;

        return copy;
    }

    public DVRP copy(int startOfNewEpoch) {
        DVRP copy = this.copy();
        copy.increaseEpoche(startOfNewEpoch);
        return copy;
    }

    public boolean isValid() {
        return vehicles.stream().allMatch(Vehicle::isValidRoute) && vehicles.size() <= noVehicles;
    }

    public Set<Integer> getUnroutedCustomers() {
        return unroutedCustomers;
    }

    public List<Vehicle> getVehicles() {
        return vehicles;
    }

    public boolean isUnusedVehicleInDepot() {
        return noVehicles > vehicles.size();
    }

    public Vehicle getNewVehicle() {
        Vehicle v = new Vehicle(this, vehicleCapacity, vehicleFixCost);
        vehicles.add(v);
        return v;
    }

    public double getFunctionValue() {
        return getFunctionValue(false);
    }

    public double getFunctionValue(boolean ingnoreUnrouted) {
        long cost = vehicles.stream().mapToLong(Vehicle::getCost).sum();

        if (vehicles.size() > noVehicles)
            cost += DataModel.lambda * getNodesInRoutes().size() * (vehicles.size() - noVehicles) * DataModel.functionPenaltyValue;

        if (unroutedCustomers.size() > 0 && !ingnoreUnrouted)
            cost += DataModel.lambda * unroutedCustomers.size() * DataModel.functionPenaltyValue;

        return cost;
    }

    // returns all nodes that can still be modified in solution
    public List<Integer> getNodesInRoutes() {
        return vehicles.stream().map(Vehicle::getRoute).flatMap(List::stream).collect(Collectors.toList());
    }

    public void increaseEpoche(int startOfNewEpoch) {
        vehicles.forEach(v -> v.increaseEpoche(startOfNewEpoch));
        addCustomerRequestInInterval(epochStartTime, startOfNewEpoch);
        epochStartTime = startOfNewEpoch;

    }

    public boolean swapCustomers(int first, int second) {
        if (first == second)
            return false;

        int swappedNodes = 0;
        Iterator<Vehicle> it = getVehicles().iterator();
        while (swappedNodes < 2) {
            Vehicle v = it.next();
            swappedNodes += v.swapCustomers(first, second);
        }

        return isValid();
    }

    public void removeCustomerFromRoute(int node) {
        Iterator<Vehicle> it = getVehicles().iterator();
        Vehicle v = it.next();
        while (!v.removeFromRoute(node)) {
            v = it.next();
        }

        if (v.getRoute().isEmpty() && v.getCompletedRoute().stream().noneMatch(x -> x > 0))
            getVehicles().remove(v);

        unroutedCustomers.add(node);
    }

    public void removeCustomersFromRoute(List<Integer> node) {
        for (Integer n : node) {
            removeCustomerFromRoute(n);
        }
    }

    public boolean newRouteForCustomer(int node) {
        removeCustomerFromRoute(node);
        Vehicle v = getNewVehicle();
        v.addToRoute(node);
        unroutedCustomers.remove(Integer.valueOf(node));
        return isValid();
    }

    public Vehicle getVehicle(int index) {
        return vehicles.get(index);
    }

    public int getCustomerCount() {
        return vehicles.stream().mapToInt(v -> v.getCompletedRoute().size() + v.getRoute().size()).sum();
    }

    public boolean removeUnroutedCustomer(int customer) {
        return unroutedCustomers.remove(customer);
    }

    public void addCustomerBeforeOther(int nodeToInsert, int insertBeforeNode) {
        removeCustomerFromRoute(nodeToInsert);
        Iterator<Vehicle> it = getVehicles().iterator();
        Vehicle v = it.next();
        while (!v.insertBeforeNode(nodeToInsert, insertBeforeNode)) {
            v = it.next();
        }
        unroutedCustomers.remove(nodeToInsert);
    }

    public List<Integer> getCompletedNodes() {
        return vehicles.stream().map(Vehicle::getCompletedRoute).flatMap(List::stream).filter(v -> v > 0).collect(Collectors.toList());
    }

    public int getTotalNodeNumber() {
        Set<Integer> s = new HashSet<>();
        s.addAll(getUnroutedCustomers());
        s.addAll(getNodesInRoutes());
        s.addAll(getCompletedNodes());
        return s.size();
    }

    public void removeUnvisited() {
        getVehicles().stream().forEach(Vehicle::clearRoute);
    }
}
