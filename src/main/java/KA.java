import static java.lang.Math.*;

/**
 * Created by Vavi on 02.04.2017.
 */
public class KA {
    public static double Rysl = 6356.766; //Условный радиус Земли
    public static double g0 = 9.80665;//ускор своб пад
    public static double R = 287.05287;//Газ пост
    public static double g0R = 0.0341632188;
    public static double Rekv = 6378.1; //Экваториальный радиус Земли
    public static double fcg = 1/298.258;//Сжатие Земли
    public static double u = 398600;//грав пост
    public static double cx = 1.3;//коэф силы лоб сопр
    public static double Smid = 15; //площадь миделя
    public static double m = 8000; // масса КА
    public static double w = 7.292115855;//Угловая скорость вращ Земли
    public static double K = 0.3; // Аэродинамическое качество
    //5 вариант
    public static double x = 4206.138;
    public static double y = 2914.179;
    public static double z = 3960.598;
    public static double Vx = -5.918745;
    public static double Vy = 3.008711;
    public static double Vz = 3.746175;
    public static double yvx = -1.5;// Угол входа
    //Таблица для плотности
    public static double[][] Plot =
            {{0,288.15,-0.0065,1.24915236*Math.pow(10,-1)},
                    {11,216.65,0,3.71093080*Math.pow(10,-1)},
                    {20,216.65,0.0010,8.97702069 *Math.pow(10,-3)},
                    {32,228.65,0.0028,1.34856449*Math.pow(10,-3)},
                    {47,270.65,0,1.45566528 *Math.pow(10,-4)},
                    {51,270.65,-0.0028,8.78587489 *Math.pow(10,-5)},
                    {71,214.65,-0.0020,6.54764879 *Math.pow(10,-6)},
                    {80,196.65,-0.0020,1.60099524 *Math.pow(10,-6)},
                    {85,186.65,0,6.91423676 *Math.pow(10,-7)},
                    {94,186.65,0.0030,1.33266712 *Math.pow(10,-7)},
                    {102.45,212.00,-0.0110,2.74594280 *Math.pow(10,-8)},
                    {117.77,380.60,0.0078,2.48852564 *Math.pow(10,-9)}};
    public static double p = 0.0;//Плотность нижней границы
    public static double am = 0.0;//Градниент температуры
    public static double Tm = 0.0;//абсолют температура
    public static double Fniz = 0.0;//Нижняя граница высоты
    public static void main(String args[]){
        double h=1;
        double a = 5;
        double b = 100;
        Double n = (b-a)/h;
        double[] X = new double[n.intValue()+1];
        double[] Y = new double[n.intValue()+1];
        int j=0;
//        for(int i=100;i>=5;i--){
//            X[j] = i;
//
//            Y[j] =PO(i);
//            j++;
//        }
//        for(int i=0;i<X.length;i++){
//            System.out.print(X[i]+" ");
//        }
//        System.out.print("\n");
//        for(int i=0;i<Y.length;i++){
//            System.out.print(Y[i]+" ");
//        }
        System.out.println("Hello");
    }
    public static double[] Fun(double[] vect){
        double[] Exit = new double[vect.length];
        double x = vect[0];
        double y = vect[1];
        double z = vect[2];
        double[] xyz = {x,y,z};

        double Vx = vect[3];
        double Vy = vect[4];
        double Vz = vect[5];
        double[] VV = {Vx,Vy,Vz};

        double r = Math.sqrt(Math.pow(x,2)+Math.pow(y,2)+Math.pow(z,2));
        double V = Math.sqrt(Math.pow(Vx,2)+Math.pow(Vy,2)+Math.pow(Vz,2));
        double[] VR = VectProizv(VV,xyz);
        double VRP = Math.sqrt(Math.pow(VR[0],2)+Math.pow(VR[1],2)+Math.pow(VR[2],2));
        for(int i=0; i<VR.length;i++){
            VR[i] = VR[i]/VRP;
        }
        double[] mm = VectProizv(VR,VV);
        double mmm = Math.sqrt(Math.pow(mm[0],2)+Math.pow(mm[1],2)+Math.pow(mm[2],2));
        for(int i=0; i<mm.length;i++) {
            mm[i] = mm[i] / mmm;
        }
        double m1 = mm[0];
        double m2 = mm[1];
        double m3 = mm[2];
        double l1 = VR[0];
        double l2 = VR[1];
        double l3 = VR[2];
        Exit[0] = Vx;
        Exit[1] = Vy;
        Exit[2] = Vz;
        Exit[3] = (-u/Math.pow(r,3))*x-cx*(Smid/(2*m))*PO(xyz)*V*Vx+
                K*cx*(Smid/(2*m))*PO(xyz)*Math.pow(V,2)*(m1*Math.cos(yvx*0.0174533)+l1*Math.sin(yvx*0.0174533))
                +x*Math.pow(w,2)+2*w*Vy;

        Exit[4] = (-u/Math.pow(r,3))*y-cx*(Smid/(2*m))*PO(xyz)*V*Vy+
                K*cx*(Smid/(2*m))*PO(xyz)*Math.pow(V,2)*(m2*Math.cos(yvx*0.0174533)+l2*Math.sin(yvx*0.0174533))
                +y*Math.pow(w,2)-2*w*Vx;

        Exit[5] = (-u/Math.pow(r,3))*z-cx*(Smid/(2*m))*PO(xyz)*V*Vz+
                K*cx*(Smid/(2*m))*PO(xyz)*Math.pow(V,2)*(m3*Math.cos(yvx*0.0174533)+l3*Math.sin(yvx*0.0174533));
        return Exit;
    }
    public static double PO(double[] xyz){//double[] xyz
        double x = xyz[0];
        double y= xyz[1];
        double z = xyz[2];
        double h = Math.sqrt(Math.pow(x,2)+Math.pow(y,2)+Math.pow(z,2))-
                Rekv*(1-(fcg*Math.pow(z,2))/(Math.pow(x,2)+Math.pow(y,2)+Math.pow(z,2)));
        double F = h*Rysl/(h+Rysl);

        Tab(F);
        double T = Tm + am*(F-Fniz);
        double po = 0.0;
        if(am!=0){
            po = p*Math.exp(1+g0R*(1/am)*Math.log(Tm/T));
        }else{
            po = p*Math.exp(g0R*(Fniz*F/Tm));
        }
        return po;
    }
    public static double Plotnos(){

        return 0.0;
    }
    public static void Tab(double F){
        for(int i=0;i<Plot.length-1;i++){
            if((F>=Plot[i][0])&(F<=Plot[i+1][0])){
                Fniz = Plot[i][0];
                p = Plot[i][0]*g0;//Плотность нижней границы
                am = Plot[i][2];//Градниент температуры
                Tm = Plot[i][1];//абсолют температура
            }
        }
    }
    public static double[] VectProizv(double[] A, double[] B){
        double[] Exit = new double [A.length];
        for(int i=0;i<3;i++){
            Exit[i] = A[(i + 1) % 3] * B[(i + 2) % 3] - A[(i + 2) % 3] * B[(i + 1) % 3];
        }
        return Exit;
    }
}
