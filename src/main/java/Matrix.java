import javax.swing.*;

/**
 * Created by Vavi on 02.04.2017.
 */
public class Matrix {
    //Нахождение обратной матрицы
    public static final void invert(double A[][]) {
        int n = A.length;
        int row[] = new int[n];
        int col[] = new int[n];
        double temp[] = new double[n];
        int hold, I_pivot, J_pivot;
        double pivot, abs_pivot;

        if (A[0].length != n) {
            System.out.println("Error in Matrix.invert, inconsistent array sizes.");
        }
        // установиим row и column как вектор изменений.
        for (int k = 0; k < n; k++) {
            row[k] = k;
            col[k] = k;
        }
        // начало главного цикла
        for (int k = 0; k < n; k++) {
            // найдем наибольший элемент для основы
            pivot = A[row[k]][col[k]];
            I_pivot = k;
            J_pivot = k;
            for (int i = k; i < n; i++) {
                for (int j = k; j < n; j++) {
                    abs_pivot = Math.abs(pivot);
                    if (Math.abs(A[row[i]][col[j]]) > abs_pivot) {
                        I_pivot = i;
                        J_pivot = j;
                        pivot = A[row[i]][col[j]];
                    }
                }
            }
            if (Math.abs(pivot) < 1.0E-10) {
                System.out.println("Matrix is singular !");
                return;
            }
            //перестановка к-ой строки и к-ого столбца с стобцом и строкой, содержащий основной элемент(pivot основу)
            hold = row[k];
            row[k] = row[I_pivot];
            row[I_pivot] = hold;
            hold = col[k];
            col[k] = col[J_pivot];
            col[J_pivot] = hold;
            // k-ую строку с учетом перестановок делим на основной элемент
            A[row[k]][col[k]] = 1.0 / pivot;
            for (int j = 0; j < n; j++) {
                if (j != k) {
                    A[row[k]][col[j]] = A[row[k]][col[j]] * A[row[k]][col[k]];
                }
            }
            // внутренний цикл
            for (int i = 0; i < n; i++) {
                if (k != i) {
                    for (int j = 0; j < n; j++) {
                        if (k != j) {
                            A[row[i]][col[j]] = A[row[i]][col[j]] - A[row[i]][col[k]] *
                                    A[row[k]][col[j]];
                        }
                    }
                    A[row[i]][col[k]] = -A[row[i]][col[k]] * A[row[k]][col[k]];
                }
            }
        }
        // конец главного цикла

        // переставляем назад rows
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                temp[col[i]] = A[row[i]][j];
            }
            for (int i = 0; i < n; i++) {
                A[i][j] = temp[i];
            }
        }
        // переставляем назад columns
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                temp[row[j]] = A[i][col[j]];
            }
            for (int j = 0; j < n; j++) {
                A[i][j] = temp[j];
            }
        }
    }
    //Умножение матрицы на вектор
    public static final void multiply(double A[][], double B[], double C[]) {
        int n = A.length;
        for (int i = 0; i < n; i++) {
            C[i] = 0.0;
            for (int j = 0; j < n; j++) {
                C[i] = C[i] + A[i][j] * B[j];
            }
        }
    }
    //Разность векторов
    public static final void raz(double A[], double B[], double C[]) {
        for (int i = 0; i < 3; i++) {
            C[i] = A[i]-B[i];

        }
    }
    //Векторное произведение
    public static double VectorProizv(double [] arr1, double [] arr2) {
        double rez = 0;
        for (int i=0;i<arr1.length;i++){
            rez = rez + arr1[i]*arr2[i];
        }
        return rez;
    }
    //Произведение двух матриц
    public static double[][] MatrPro(double[][] A,double[][] B){
        double[][] Exit = null;
        int m = A.length;
        int n = B[0].length;
        int o = B.length;
        if(A[0].length!=B.length){
            JOptionPane.showMessageDialog(null, "Неправильные матрицы!");
        }else{
            Exit = new double[m][n];
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    for (int k = 0; k < o; k++) {
                        Exit[i][j] += A[i][k] * B[k][j];
                    }
                }
            }
        }
        return Exit;
    }
    //Транспонирование  матрицы
    public static double[][] transpose(double[][] matrix) throws IllegalArgumentException {
        if(matrix.length == 0) {
            throw new IllegalArgumentException("Empty array");
        }
        int rowLength = matrix[0].length;
        for (double[] ai:matrix) {
            if (rowLength != ai.length) {
                throw new IllegalArgumentException("Non-equal rows");
            }
        }

        double [][] tMatrix = new double[rowLength][];
        for (int i = 0; i < rowLength; i++) {
            tMatrix[i] = new double[matrix.length];
        }
        for (int i = 0; i < matrix.length; i++) {
            double[] tArr = matrix[i];
            for (int j = 0; j < rowLength; j++) {
                tMatrix[j][i] = tArr[j];
            }
        }
        return tMatrix;
    }
    //Произведение кватернионов 4 порядка
    public static final void proizkvat(double B[], double A[], double C[]) {

        C[0] = A[0]*B[0] - A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
        C[1] = A[0]*B[1] + A[1]*B[0] + A[3]*B[2] - A[2]*B[3];
        C[2] = A[0]*B[2] + A[2]*B[0] + A[1]*B[3] - A[3]*B[1];
        C[3] = A[0]*B[3] + A[3]*B[0] + A[2]*B[1] - A[1]*B[2];
    }
    //Умножение вектора на число
    public static final void vekchi(double A[], double B, double C[]) {
        for (int i = 0; i < A.length; i++) {
            C[i] = A[i]*B;
        }
    }
    //Нормировка вектора
    public static final void norm(double q[]) {
        double E = 1.0/(Math.sqrt(Matrix.VectorProizv(q,q)));
        Matrix.vekchi(q, E, q);

    }
    //Добавление нулевого столбца
    public static final void zero(double q[],double q2[]) {
        q2[0]=0;
        for(int i=0;i<3;i++){
            q2[i+1]=q[i];
        }

    }
    //Линеарезация
    public static double linrez(double x1,double y1, double x2,double y2,double x) {
        double y=y1+((y2-y1)*(x-x1))/(x2-x1);
        return y;
    }
}