package sczuka.tech;

import java.util.ArrayList;
import java.util.List;

class Vehicle {
    private DVRP dvrp;
    private List<Integer> route = new ArrayList<>();
    private List<Integer> completedRoute = new ArrayList<>(List.of(0));
    private int timeOfCompletedRoute = 0;
    private int noTimeViolationsCompletedRoute = 0;
    private int currentDemand = 0;
    private int currentTime = 0;
    private int capacity = 0;
    private int fixCost = 0;
    private int noTimeViolations = 0;

    Vehicle(DVRP dvrp, int capacity, int fixCost) {
        this.dvrp = dvrp;
        this.capacity = capacity;
        this.fixCost = fixCost;
    }

    public Vehicle copy(DVRP dvrp) {
        Vehicle copy = new Vehicle(dvrp, capacity, fixCost);
        copy.route = new ArrayList<>(route);
        copy.completedRoute = new ArrayList<>(completedRoute);
        copy.timeOfCompletedRoute = timeOfCompletedRoute;
        copy.noTimeViolationsCompletedRoute = noTimeViolationsCompletedRoute;
        copy.currentDemand = currentDemand;
        copy.currentTime = currentTime;
        copy.noTimeViolations = noTimeViolations;

        return copy;
    }

    public int getPosition() {
        if (route.isEmpty())
            return completedRoute.get(completedRoute.size() - 1);
        else
            return route.get(route.size() - 1);
    }

    public List<Integer> getRoute() {
        return route;
    }

    public List<Integer> getCompletedRoute() {
        return completedRoute;
    }

    public int getStartPostion() {
        return completedRoute.get(completedRoute.size() - 1);
    }

    public void increaseEpoche(int epochStartTime) {
        if (route.isEmpty())
            return;

        while (timeOfCompletedRoute <= epochStartTime && !route.isEmpty()) {
            timeOfCompletedRoute += DataModel.timeMatrix[getStartPostion()][route.get(0)];
            if (timeOfCompletedRoute > DataModel.timeWindows[route.get(0)][1]) {
                noTimeViolationsCompletedRoute++;
            }
            completedRoute.add(route.remove(0));
        }
    }

    public boolean isValidRoute(int j) {
        return currentDemand + DataModel.demands[j] <= capacity && noTimeViolations == 0
                && currentTime + DataModel.timeMatrix[getPosition()][j] <= DataModel.timeWindows[j][1];
    }

    public boolean isValidRoute() {
        return currentDemand <= capacity && noTimeViolations == 0 && noTimeViolationsCompletedRoute == 0;
    }

    public boolean addToRoute(int j) {
        if (dvrp != null)
            dvrp.removeUnroutedCustomer(j);

        currentDemand += DataModel.demands[j];
        currentTime = Math.max(currentTime + DataModel.timeMatrix[getPosition()][j], DataModel.timeWindows[j][0]);
        if (currentTime > DataModel.timeWindows[j][1]) {
            noTimeViolations++;
        }

        route.add(j);
        return isValidRoute();
    }

    public boolean addToRoute(List<Integer> list) {
        for (int i : list) {
            addToRoute(i);
        }
        return isValidRoute();
    }

    public boolean addToRoute(int position, int j) {
        if (dvrp != null)
            dvrp.removeUnroutedCustomer(j);

        currentDemand += DataModel.demands[j];
        route.add(position, j);
        updateTime();
        return isValidRoute();
    }

    public boolean addToRoute(int position, List<Integer> list) {
        if (dvrp != null)
            dvrp.getUnroutedCustomers().removeAll(list);

        List<Integer> newRoute = new ArrayList<>();
        if (position > 0)
            newRoute.addAll(getRoute().subList(0, position));

        newRoute.addAll(list);

        if (position < getRoute().size())
            newRoute.addAll(getRoute().subList(position, getRoute().size()));

        route = newRoute;

        currentDemand += list.stream().mapToInt(j -> DataModel.demands[j]).sum();
        updateTime();
        return isValidRoute();
    }

    public boolean removeFromRoute(int j) {
        if (route.remove(Integer.valueOf(j))) {
            currentDemand -= DataModel.demands[j];
            updateTime();
            return true;
        } else {
            return false;
        }
    }

    private void updateTime() {
        currentTime = timeOfCompletedRoute;
        noTimeViolations = 0;

        if (!route.isEmpty()) {
            currentTime += DataModel.timeMatrix[getStartPostion()][route.get(0)];
            if (currentTime > DataModel.timeWindows[route.get(0)][1]) {
                noTimeViolations++;
            }

            for (int j = 0; j < route.size() - 1; j++) {
                currentTime = Math.max(currentTime + DataModel.timeMatrix[route.get(j)][route.get(j + 1)], DataModel.timeWindows[route.get(j + 1)][0]);
                if (currentTime > DataModel.timeWindows[route.get(j + 1)][1]) {
                    noTimeViolations++;
                }
            }
        }
    }

    public int swapCustomers(int i, int j) {
        Integer first = Integer.valueOf(i);
        Integer second = Integer.valueOf(j);

        int indexOfFirst = route.indexOf(first);
        int indexOfSecond = route.indexOf(second);

        int replacedValues = 0;
        if (indexOfFirst != -1) {
            route.set(indexOfFirst, second);
            currentDemand -= DataModel.demands[first];
            currentDemand += DataModel.demands[second];
            replacedValues++;
        }
        if (indexOfSecond != -1) {
            route.set(indexOfSecond, first);
            currentDemand += DataModel.demands[first];
            currentDemand -= DataModel.demands[second];
            replacedValues++;
        }

        if (replacedValues > 0)
            updateTime();

        return replacedValues;

    }

    public boolean insertBeforeNode(int nodeToInsert, int insertBeforeNode) {
        int index = route.indexOf(insertBeforeNode);
        if (index == -1)
            return false;

        route.add(index, nodeToInsert);
        currentDemand += DataModel.demands[nodeToInsert];
        updateTime();
        return true;

    }

    public long getCost() {
        if (route.isEmpty())
            return fixCost + timeOfCompletedRoute;

        long cost = fixCost + timeOfCompletedRoute + DataModel.timeMatrix[getStartPostion()][route.get(0)]
                + DataModel.timeMatrix[route.get(route.size() - 1)][0];
        for (int i = 0; i < route.size() - 1; i++) {
            cost += DataModel.timeMatrix[route.get(i)][route.get(i + 1)];
        }

        if (!isValidRoute())
            cost += DataModel.lambda * getNoViolations() * DataModel.functionPenaltyValue;

        return cost;
    }

    public int getNoViolations() {
        return (currentDemand > capacity ? 1 : 0) + noTimeViolations + noTimeViolationsCompletedRoute;
    }

    public void clearRoute() {
        dvrp.getUnroutedCustomers().addAll(getRoute());
        List<Integer> removed = new ArrayList<>(getRoute());
        removed.forEach(r -> {
            currentDemand -= DataModel.demands[r];
        });
        getRoute().clear();
        updateTime();
    }

    public List<Integer> removeRouteBetweenIndex(int i1, int i2) {
        List<Integer> temp = new ArrayList<>(getRoute());

        List<Integer> removedCustomers = new ArrayList<>();

        if (i1 < i2)
            for (int i = i1; i < i2 + 1; i++) {
                removedCustomers.add(temp.get(i));
                removeFromRoute(temp.get(i));
            }
        else
            for (int i = i2; i < i1 + 1; i++) {
                removedCustomers.add(temp.get(i));
                removeFromRoute(temp.get(i));
            }
        return removedCustomers;
    }
}
