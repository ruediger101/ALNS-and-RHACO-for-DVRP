package sczuka.tech;

public class TriangleCalculations {

    public static double calculateArea(double a, double b, double c) {
        // Using Heron's formula to calculate area
        double s = (a + b + c) / 2;
        return Math.sqrt(s * (s - a) * (s - b) * (s - c));
    }

    public static double calculateHeight(double area, double base) {
        // Calculate height from area and base length
        return (2 * area) / base;
    }

    public static double[] calculateAngles(double a, double b, double c) {
        // Calculate angles using the law of cosines
        double angleA = Math.acos((b * b + c * c - a * a) / (2 * b * c));
        double angleB = Math.acos((a * a + c * c - b * b) / (2 * a * c));
        double angleC = Math.acos((a * a + b * b - c * c) / (2 * a * b));

        // Converting angles from radians to degrees
        return new double[]{Math.toDegrees(angleA), Math.toDegrees(angleB), Math.toDegrees(angleC)};
    }

    public static void main(String[] args) {
        // Example distances between points A, B, and C
        double a = 2.12836; // Distance BC
        double b = 4; // Distance AC
        double c = 5.38919; // Distance AB (Hypotenuse in this right triangle example)

        double area = calculateArea(a, b, c);
        double height = calculateHeight(area, a);
        double[] angles = calculateAngles(a, b, c);

        System.out.println("Area of the triangle: " + area);
        System.out.println("Height of the triangle relative to side BC: " + height);
        System.out.println("Angle at A: " + angles[0] + " degrees");
        System.out.println("Angle at B: " + angles[1] + " degrees");
        System.out.println("Angle at C: " + angles[2] + " degrees");
    }
}


