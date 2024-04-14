package sczuka.tech;

public class ClosestPointOnSegment {
    static class Point {
        double x;
        double y;

        public Point(double x, double y) {
            this.x = x;
            this.y = y;
        }
        
        // Method to subtract two points
        public Point subtract(Point p) {
            return new Point(x - p.x, y - p.y);
        }
        
        // Method to add two points
        public Point add(Point p) {
            return new Point(x + p.x, y + p.y);
        }
        
        // Method to multiply point by a scalar
        public Point multiply(double scalar) {
            return new Point(x * scalar, y * scalar);
        }
        
        // Dot product of two points
        public double dot(Point p) {
            return x * p.x + y * p.y;
        }
        
        // Compute the distance to another point
        public double distanceTo(Point p) {
            return Math.sqrt((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y));
        }
    }
    
    public static Point closestPoint(Point A, Point B, Point P) {
        Point AB = B.subtract(A);
        Point AP = P.subtract(A);
        
        // Project point P onto the line segment
        double ab2 = AB.dot(AB);
        double ap_ab = AP.dot(AB);
        double t = ap_ab / ab2;

        // Clamping t to the range [0,1] will ensure the point lies on the segment
        t = Math.max(0, Math.min(1, t));

        // Calculate the closest point
        Point closest = A.add(AB.multiply(t));
        return closest;
    }
    
    public static void main(String[] args) {
        Point A = new Point(0, 3);
        Point B = new Point(0, 4);
        Point P = new Point(3, 2);

        Point closest = closestPoint(A, B, P);
        System.out.println("The closest point is: (" + closest.x + ", " + closest.y + ")");
    }
}
