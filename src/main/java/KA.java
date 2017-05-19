import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import static java.lang.Math.*;

/**
 * Created by Vavi on 02.04.2017.
 */
public class KA {
    public static double Rysl = 6356.766; //Условный радиус Земли метры км
    public static double RZ = 6371; //радиус Земли метры км
    public static double g0 = 9.80665 *Math.pow(10,-3);//ускор своб пад км
    public static double R = 287.05287*Math.pow(10,-3);//Газ пост
    public static double g0R = (g0/R)*Math.pow(10,3);//км
    public static double Rekv = 6378.1; //Экваториальный радиус Земли км
    public static double fcg = 0.003352824419;//Сжатие Земли
    public static double u = 398600;//грав пост км
    public static double cx = 1.3;//коэф силы лоб сопр
    public static double Smid = 15 / Math.pow(1000, 2); //площадь миделя км
    public static double m = 8000; // масса КА
    public static double w = 7.292115855 * Math.pow(10, -5);//Угловая скорость вращ Земли
    public static double K = 0.3; // Аэродинамическое качество
    //5 вариант
    public static double x = 4206.138;
    public static double y = 2914.179;
    public static double z = 3960.598;
    public static double Vx = -5.918745;
    public static double Vy = 3.008711;
    public static double Vz = 3.746175;
    //6 вариант
//    public static double x = 3701.992;
//    public static double y = 3158.130;
//    public static double z = 4265.507;
//    public static double Vx = -6.380016;
//    public static double Vy = 2.661533;
//    public static double Vz = 3.164187;
    //Вариант хз
//    public static double x = 5065.124;
//    public static double y = 2373.486;
//    public static double z = 3253.021;
//    public static double Vx = -4.857479;
//    public static double Vy = 3.565160;
//    public static double Vz = 4.695924;

    public static double yvx = 45;// Угол входа
    //Таблица для плотности
    public static double[][] Plot =
            {{0,     288.15,  -0.0065,  1.24915236 *Math.pow(10,-1)},
                    {11,     216.65,  -0.0065,  3.71093080 *Math.pow(10,-2)},
                    {20,     216.65,   0,       8.97702069 *Math.pow(10,-3)},
                    {32,     228.65,   0.0010,  1.34856449 *Math.pow(10,-3)},
                    {47,     270.65,   0.0028,  1.45566528 *Math.pow(10,-4)},
                    {51,     270.65,   0,       8.78587489 *Math.pow(10,-5)},
                    {71,     214.65,  -0.0028,  6.54764879 *Math.pow(10,-6)},
                    {80,     196.65,  -0.0017,  1.60099524 *Math.pow(10,-6)},
                    {85,     186.65,  -0.0020,  6.91423676 *Math.pow(10,-7)},
                    {94,     186.65,   0,       1.33266712 *Math.pow(10,-7)},
                    {102.45, 212.00,   0.0030,  2.74594280 *Math.pow(10,-8)},
                    {117.77, 380.60,   0.0110,  2.48852564 *Math.pow(10,-9)}};
    public static double p = 0.0;//Плотность нижней границы
    public static double am = 0.0;//Градниент температуры
    public static double Tm = 0.0;//абсолют температура
    public static double Fniz = 0.0;//Нижняя граница высоты
    public static double F = 0.0;//Высота КА
    public static double Fgeom = 0.0;//Высота КА
    public static double P = 0.0;//Давление атмосферы

    public static void main(String args[]) {
        for (int i = 0; i < Plot.length; i++) {
            Plot[i][2] = Plot[i][2] *1000;
            Plot[i][3] = Plot[i][3] * Math.pow(10, 12);
        }
        double h = 0.2;
        double time0 = 0;
        double h1 = Math.sqrt(Math.pow(x, 2) + Math.pow(y, 2) + Math.pow(z, 2)) -
                Rekv * (1 - (fcg * Math.pow(z, 2)) / (Math.pow(x, 2) + Math.pow(y, 2) + Math.pow(z, 2)));
        double visota0 = h1;

        ArrayList<Double> time = new ArrayList<Double>();
        ArrayList<Double> Visota = new ArrayList<Double>();
        ArrayList<Double> Atmor = new ArrayList<Double>();
        ArrayList<Double> LX = new ArrayList<Double>();
        ArrayList<Double> LY = new ArrayList<Double>();
        ArrayList<Double> LZ = new ArrayList<Double>();
        ArrayList<Double> LVX = new ArrayList<Double>();
        ArrayList<Double> LVY = new ArrayList<Double>();
        ArrayList<Double> LVZ = new ArrayList<Double>();
        ArrayList<Double> LV = new ArrayList<Double>();
        ArrayList<Double> n = new ArrayList<Double>();
        ArrayList<double[]> GeogrKor = new ArrayList<double[]>();
        ArrayList<double[]> ProvSost = new ArrayList<double[]>();
        ArrayList<double[]> ProvBuf = new ArrayList<double[]>();
        ArrayList<Integer> ProvSostV = new ArrayList<Integer>();
        ArrayList<Double> ProstoP = new ArrayList<Double>();
        time.add(time0);
        Visota.add(visota0);
        Atmor.add(0.0);
        n.add(0.0);
        int i = 0;
        double[] vect = {x, y, z, Vx, Vy, Vz};
        double[] Fil = GetFila(new double[]{vect[0],vect[1],vect[2]});
        GeogrKor.add(new double[]{Fil[0]*(180/Math.PI),Fil[1]*(180/Math.PI)});
        LX.add(x);
        LY.add(y);
        LZ.add(z);
        LVX.add(Vx);
        LVY.add(Vy);
        LVZ.add(Vz);
        double v = Math.sqrt(Math.pow(Vx, 2) + Math.pow(Vy, 2) + Math.pow(Vz, 2));
        LV.add(v);
//        double k = 100;
//        while(k>0){
//            double p = Fpo(k);
//            Visota.add(k);
//            Atmor.add(p);
//            k=k-0.035;
//        }
//        String str = "";
//        for(int ii=0;ii<Visota.size();ii++){
//            str = str +Visota.get(ii)+";"+Atmor.get(ii)+"\n";
//        }
//        TextIn(str,"testing.txt");
        int o = 50;
        double[] buf = null;
        double[] tmp = null;
        while (visota0 > 5) {
            buf = vect;
            //vect = SumVec(vect, UmChis(Fun(vect), h));
            vect = runge(vect,45);
            if((visota0<o)&(o>=10)){
                ProvSost.add(vect);
                ProvBuf.add(buf);
                ProvSostV.add(o);
//                System.out.println("Высота = "+o);
//                System.out.println("Высота реал = "+visota0);
                // System.out.println(vect[0]+" "+vect[1]+" "+vect[2]);
                o = o - 5;
            }
            tmp = vect;

            Fil = GetFila(new double[]{vect[0],vect[1],vect[2]});
            GeogrKor.add(new double[]{Fil[0]*(180/Math.PI),Fil[1]*(180/Math.PI)});
            if(visota0<20){
                //System.out.println(Fil[0]*(180/Math.PI)+";"+Fil[1]*(180/Math.PI));
            }
            visota0 = Fgeom;
            Visota.add(visota0);
            LX.add(vect[0]);
            LY.add(vect[1]);
            LZ.add(vect[2]);
            LVX.add(vect[3]);
            LVY.add(vect[4]);
            LVZ.add(vect[5]);
            v = Math.sqrt(Math.pow(vect[3], 2) + Math.pow(vect[4], 2) + Math.pow(vect[5], 2));
            LV.add(v);
            //System.out.println("Visota = "+visota0);
            time0 = time0 + 0.2;
            time.add(time0);
            P = PO(new double[]{vect[0],vect[1],vect[2]});
            Atmor.add(P);
            i++;
            double N = cx*Math.sqrt(1+Math.pow(K,2))*((P*Math.pow(v,2)*Smid)/(2*m*g0));
            n.add(N);
//            if(i==145){
//                System.out.println();
//            }
        }
        String str = "";
        for (int ii = 0; ii < Visota.size(); ii++) {
            str = str + time.get(ii) + ";" + Visota.get(ii) + ";" + Atmor.get(ii) + ";" + LX.get(ii)+ ";" + LY.get(ii)+ ";" + LZ.get(ii) + ";"+ LV.get(ii) + ";" + n.get(ii)+ "\n";
        }
        TextIn(str, "testing.txt");
//        ArrayList<double[]> GeogrKor = new ArrayList<double[]>();
//        for(double[] ve:ProvSost){
//            double X = ve[0];
//            double Y = ve[1];
//            double Z = ve[2];//* (Math.PI / 180)
//            double A = Math.pow((1-fcg),2);
//            double B = Math.sqrt(Math.pow(X,2)+Math.pow(Y,2));
//            double fi = Math.asin(Z/(Math.sqrt(Math.pow(A*B,2)+Math.pow(Z,2))));
//            double lam = Math.atan(Y/X);
//            GeogrKor.add(new double[]{fi*(180/Math.PI),lam*(180/Math.PI)});
//        }
        str = "";
        for (int ii = 0; ii < GeogrKor.size(); ii++) {
            str = str +GeogrKor.get(ii)[0]+";"+GeogrKor.get(ii)[1]+ "\n";
        }
        TextIn(str, "testing2.txt");
        //Строим области
        double[] Ugl = new double[]{-180, -150, -120, -90, -60, -30, -15, 0,15, 30, 60, 90, 120, 150, 180};
        str = "";
        int ii=0;
        //double[] ve_ = tmp;
        //double[] Dima = new double[]{1644.3154160063766, 3841.0247677715925, 4868.772937918053, -4.2695073602414535, 1.2171817484444023, 0.30248489408059825};
        //ФИ И ЛЯМБДА от последних координат номинала
        double[] FiLa1 =new double[]{GeogrKor.get(GeogrKor.size()-1)[0]*(Math.PI/180),GeogrKor.get(GeogrKor.size()-1)[1]*(Math.PI/180)};
        FiLa1 = new double[]{FiLa1[0]*180/Math.PI,FiLa1[1]*180/Math.PI};
        //Цикл по векторам состояния высот
        //double[] Dima10 = new double[]{1571.346241894644, 3887.1571224570048, 4803.5059919596033, -0.048859589280481987, -0.063638928291775615, -0.09562669122723165};
        //System.out.println("============================================= ");
        for(double[] ve:ProvSost){
            // ve = Dima10;
            //Массив с высотами (50,45...)
            str = str+"Высота = "+ProvSostV.get(ii)+"\n";
            str = str + "Угол;Вост;Сев;n"+"\n";
            //Получение высоты

            //Цикл по углам
            for(double u:Ugl){
                //Интегрирование с другим углом
                double he = ProvSostV.get(ii);
                double[] ve_0_ = runge(ve,u);
                //double[] ve2 = ve_0_;
                double N2 = 0;
               // System.out.println("he = "+he);
                while(he >5){
                    //Расчёт высоты
                    //Интегрирование
                    ve_0_ = runge(ve_0_,u);
                    double P2 = PO(new double[]{ve_0_[0],ve_0_[1],ve_0_[2]});
                    //System.out.println(P2);
                    v =  Math.sqrt(Math.pow(ve_0_[3], 2) + Math.pow(ve_0_[4], 2) + Math.pow(ve_0_[5], 2));
                    //System.out.println(v);
                    double N = cx*Math.sqrt(1+Math.pow(K,2))*((P2*Math.pow(v,2)*Smid)/(2*m*g0));
                    if(he==ProvSostV.get(ii)){
                        N2 = N;
                    }else{
                        if(N>N2){
                            N2 = N;
                        }
                    }
                    //System.out.println("N2 = "+N2);
                    he =  Math.sqrt(Math.pow(ve_0_[0], 2) + Math.pow(ve_0_[1], 2) + Math.pow(ve_0_[2], 2)) -
                            Rekv * (1 - (fcg * Math.pow(ve_0_[2], 2)) / (Math.pow(ve_0_[0], 2) + Math.pow(ve_0_[1], 2) + Math.pow(ve_0_[2], 2)));

                }
                //Фи и лямбда полученного последнего вектора состояния
                double[] FiLa_01 = GetFila(new double[]{ve_0_[0],ve_0_[1],ve_0_[2]});//*(180/Math.PI)
                FiLa_01 = new double[]{FiLa_01[0]*180/Math.PI,FiLa_01[1]*180/Math.PI};
                //Отклонения
                double cev_0 = ((FiLa_01[0] - FiLa1[0]))*RZ;
                double vost_0 = ((FiLa_01[1]- FiLa1[1]))*RZ*Math.cos(FiLa1[0]);//*Math.PI/180)*180/Math.PI;

                //double g1 = 9.80665;
                //double P2 = PO(new double[]{ve_0_[0],ve_0_[1],ve_0_[2]});
//                System.out.println("ve_0_[0] = "+ve_0_[0]+" ve_0_[1 = "+ve_0_[1]+" ve_0_[2] = "+ve_0_[2]+"\n");

//                System.out.println("p = "+P);
                //double N = cx*Math.sqrt(1+Math.pow(K,2))*((P2*Math.pow(v,2)*Smid)/(2*m*g0));
                str = str + u+";"+vost_0/100+";"+cev_0/100+";"+N2+"\n";
            }

//            for(double u:Ugl){
//                double[] ve_0_ = runge(ve,u);
//                double[] FiLa_01 = GetFila(new double[]{ve_0_[0],ve_0_[1],ve_0_[2]});//*(180/Math.PI)
//                double cev_0 = ((FiLa_01[0] - FiLa1[0])*180/Math.PI)*RZ;
//                double vost_0 = ((FiLa_01[1]- FiLa1[1])*180/Math.PI)*RZ*Math.cos(FiLa1[0]);
//                //System.out.println("cev = "+cev_0+" vost_0 = "+vost_0+"\n");
//
//                P = PO(new double[]{ve_0_[0],ve_0_[1],ve_0_[2]});
//                double N = cx*Math.sqrt(1+Math.pow(K,2))*((P*Math.pow(v,2)*Smid)/(2*m*g0));
//                str = str + u+";"+vost_0+";"+cev_0+";"+N+"\n";//+ve2[0]+";"+ve2[1]+";"+ve2[2]+";"+ve2[3]+";"+ve2[4]+";"+ve2[5]+"\n";
////
////                //double[] ve1 =runge(vect,u);;//SumVec(ProvBuf.get(ii), UmChis(FunUgl(ProvBuf.get(ii),u), h));;// FunUgl(ve,u);
////                //double[] ve2 =SumVec(ve1, UmChis(FunUgl(ve1,u), h));;// FunUgl(ve,u);
////                double[] ve2 = runge(vect,u);
////                //ve2 = new double[]{1644.3153409213035, 3841.024535482519, 4868.772812626712, -4.270257643742819, 1.2148600947736947, 0.3012326617518336};
////               // System.out.println(ve2[0]+" "+ve2[1]+" "+ve2[2]+" "+ve2[3]+" "+ve2[4]+" "+ve2[5]+"\n");
////                double[] FiLa = GetFila(new double[]{ve2[0],ve2[1],ve2[2]});
////
////                double cev = (FiLa[0]*(180/Math.PI) - FiLa1[0]*(180/Math.PI))*RZ;
////                double vost = (  FiLa[1]*(180/Math.PI)- FiLa1[1]*(180/Math.PI))*RZ*Math.cos(FiLa1[0]);
////                P = PO(new double[]{ve2[0],ve2[1],ve2[2]});
////                //System.out.println("P = "+P);
////                v = Math.sqrt(Math.pow(ve2[3], 2) + Math.pow(ve2[4], 2) + Math.pow(ve2[5], 2));
////                double N = cx*Math.sqrt(1+Math.pow(K,2))*((P*Math.pow(v,2)*Smid)/(2*m*g0));
////                str = str + u+";"+vost+";"+cev+";"+N+"\n";//+ve2[0]+";"+ve2[1]+";"+ve2[2]+";"+ve2[3]+";"+ve2[4]+";"+ve2[5]+"\n";
//            }
            ii++;
        }
        TextIn(str, "Polet.txt");
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

    public static double[] UmChis(double[] Vect, double chis) {
        double[] Exit = new double[Vect.length];
        for (int i = 0; i < Vect.length; i++) {
            Exit[i] = Vect[i] * chis;
        }
        return Exit;
    }

    public static double[] GetFila(double[] r){
        double X = r[0];
        double Y = r[1];
        double Z = r[2];//* (Math.PI / 180)
        double A = Math.pow((1-fcg),2);
        double B = Math.sqrt(Math.pow(X,2)+Math.pow(Y,2));
        double fi = Math.asin(Z/(Math.sqrt(Math.pow(A*B,2)+Math.pow(Z,2))));
        double lam = Math.atan(Y/X);
        return new double[]{fi, lam};
    }
    public static double[] SumVec(double[] A, double[] B) {
        double[] Exit = new double[A.length];
        for (int i = 0; i < A.length; i++) {
            Exit[i] = A[i] + B[i];
        }
        return Exit;
    }

    //Запись текста в файл без диалога
    public static void TextIn(String text, String pyt) {
        File file = new File(pyt);
        try {
            if (!file.exists()) {
                file.createNewFile();
            }
            PrintWriter out = new PrintWriter(file.getAbsoluteFile());
            try {
                //Записываем текст у файл
                out.print(text);
            } finally {
                //После чего мы должны закрыть файл
                //Иначе файл не запишется
                out.close();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static double[] Fun(double[] vect) {
        double[] Exit = new double[vect.length];
        double x = vect[0];
        double y = vect[1];
        double z = vect[2];
        double[] xyz = {x, y, z};

        double Vx = vect[3];
        double Vy = vect[4];
        double Vz = vect[5];
        double[] VV = {Vx, Vy, Vz};

        double r = Math.sqrt(Math.pow(x, 2) + Math.pow(y, 2) + Math.pow(z, 2));
        double V = Math.sqrt(Math.pow(Vx, 2) + Math.pow(Vy, 2) + Math.pow(Vz, 2));
        double[] VR = VectProizv(VV, xyz);
        double VRP = Math.sqrt(Math.pow(VR[0], 2) + Math.pow(VR[1], 2) + Math.pow(VR[2], 2));
        for (int i = 0; i < VR.length; i++) {
            VR[i] = VR[i] / VRP;
        }
        double[] mm = VectProizv(VR, VV);
        double mmm = Math.sqrt(Math.pow(mm[0], 2) + Math.pow(mm[1], 2) + Math.pow(mm[2], 2));
        for (int i = 0; i < mm.length; i++) {
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
        P = PO(xyz);
        //System.out.println("p = "+P);
        Exit[3] = ((-u / Math.pow(r, 3)) * x - cx * (Smid / (2 * m)) * P * V * Vx +
                K * cx * (Smid / (2 * m)) * P * Math.pow(V, 2) * (m1 * Math.cos(yvx) + l1 * Math.sin(yvx))
                + x * Math.pow(w, 2) + 2 * w * Vy);

        Exit[4] = ((-u / Math.pow(r, 3)) * y - cx * (Smid / (2 * m)) * P * V * Vy +
                K * cx * (Smid / (2 * m)) * P * Math.pow(V, 2) * (m2 * Math.cos(yvx) + l2 * Math.sin(yvx))
                + y * Math.pow(w, 2) - 2 * w * Vx);

        Exit[5] = ((-u / Math.pow(r, 3)) * z - cx * (Smid / (2 * m)) * P * V * Vz +
                K * cx * (Smid / (2 * m)) * P * Math.pow(V, 2) * (m3 * Math.cos(yvx) + l3 * Math.sin(yvx)));

        return Exit;
    }

    public static double[] FunUgl(double[] vect, double Ugl) {
        double[] Exit = new double[vect.length];
        double x = vect[0];
        double y = vect[1];
        double z = vect[2];
        double[] xyz = {x, y, z};

        double Vx = vect[3];
        double Vy = vect[4];
        double Vz = vect[5];
        double[] VV = {Vx, Vy, Vz};

        double r = Math.sqrt(Math.pow(x, 2) + Math.pow(y, 2) + Math.pow(z, 2));
        double V = Math.sqrt(Math.pow(Vx, 2) + Math.pow(Vy, 2) + Math.pow(Vz, 2));
        double[] VR = VectProizv(VV, xyz);
        double VRP = Math.sqrt(Math.pow(VR[0], 2) + Math.pow(VR[1], 2) + Math.pow(VR[2], 2));
        for (int i = 0; i < VR.length; i++) {
            VR[i] = VR[i] / VRP;
        }
        double[] mm = VectProizv(VR, VV);
        double mmm = Math.sqrt(Math.pow(mm[0], 2) + Math.pow(mm[1], 2) + Math.pow(mm[2], 2));
        for (int i = 0; i < mm.length; i++) {
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
        P = PO(xyz);

        //System.out.println("p = "+P);
        Exit[3] = ((-u / Math.pow(r, 3)) * x - cx * (Smid / (2 * m)) * P * V * Vx +
                K * cx * (Smid / (2 * m)) * P * Math.pow(V, 2) * (m1 * Math.cos(Ugl*(180/Math.PI)) + l1 * Math.sin(Ugl*(180/Math.PI)))
                + x * Math.pow(w, 2) + 2 * w * Vy);

        Exit[4] = ((-u / Math.pow(r, 3)) * y - cx * (Smid / (2 * m)) * P * V * Vy +
                K * cx * (Smid / (2 * m)) * P * Math.pow(V, 2) * (m2 * Math.cos(Ugl*(180/Math.PI)) + l2 * Math.sin(Ugl*(180/Math.PI)))
                + y * Math.pow(w, 2) - 2 * w * Vx);

        Exit[5] = ((-u / Math.pow(r, 3)) * z - cx * (Smid / (2 * m)) * P * V * Vz +
                K * cx * (Smid / (2 * m)) * P * Math.pow(V, 2) * (m3 * Math.cos(Ugl*(180/Math.PI)) + l3 * Math.sin(Ugl*(180/Math.PI))));

        return Exit;
    }


    public static double PO(double[] xyz) {//double[] xyz
        double x = xyz[0];
        double y = xyz[1];
        double z = xyz[2];
        Fgeom = Math.sqrt(Math.pow(x, 2) + Math.pow(y, 2) + Math.pow(z, 2)) -
                Rekv * (1 - (fcg * Math.pow(z, 2)) / (Math.pow(x, 2) + Math.pow(y, 2) + Math.pow(z, 2)));
        //System.out.println("f = "+Fgeom);
        F = Fgeom * Rysl / (Fgeom + Rysl);
        return Fpo(F);
    }

    public static double Fpo(double F) {
        Tab(F);
        double T = Tm + am * (F - Fniz);
        double po = 0.0;
        if (am != 0) {
            po = p * Math.exp((1 + g0R / am) * Math.log(Tm / T));
        } else {
            po = p * Math.exp(g0R * (Fniz - F) / Tm);
        }
        //po = po/Math.pow(1000,2);
        if (po == 0.0) {
            System.out.println();
        }
        return po;
    }

    public static void Tab(double F) {
        for (int i = 0; i < Plot.length - 1; i++) {
            if ((F >= Plot[i][0]) & (F < Plot[i + 1][0])) {
                Fniz = Plot[i+1][0];
                Tm = Plot[i+1][1];//абсолют температура
                am = Plot[i+1][2];//Градниент температуры
                p = Plot[i+1][3]*g0;//Плотность нижней границы
                break;
            }
        }
    }

    public static double[] right(double[] vect, double gam){
        double c_xa = 1.3;
        double S_mid = 15 * Math.pow(10,-6);
        double m = 8000;
        double mu = 398600;
        double omega = 7.292115855 * Math.pow(10,-5);
        double K = 0.3;

        gam = gam*(Math.PI/180);

        double[] res = new double[6];

        double r = Math.sqrt(Math.pow(vect[0], 2) + Math.pow(vect[1], 2) + Math.pow(vect[2], 2));
        double V = Math.sqrt(Math.pow(vect[3], 2) + Math.pow(vect[4], 2) + Math.pow(vect[5], 2));
        double[] xyz = {vect[0], vect[1], vect[2]};
        double[] VV =  {vect[3], vect[4], vect[5]};
        double[] VR = VectProizv(VV, xyz);
        double VRP = Math.sqrt(Math.pow(VR[0], 2) + Math.pow(VR[1], 2) + Math.pow(VR[2], 2));
        for (int i = 0; i < VR.length; i++) {
            VR[i] = VR[i] / VRP;
        }
        double[] mm = VectProizv(VR, VV);
        double mmm = Math.sqrt(Math.pow(mm[0], 2) + Math.pow(mm[1], 2) + Math.pow(mm[2], 2));
        for (int i = 0; i < mm.length; i++) {
            mm[i] = mm[i] / mmm;
        }
        double m1 = mm[0];
        double m2 = mm[1];
        double m3 = mm[2];
        double l1 = VR[0];
        double l2 = VR[1];
        double l3 = VR[2];
        double ro = PO(xyz);

        //# первые три уравнения
        res[0] = VV[0];
        res[1] = VV[1];
        res[2] = VV[2];
        double x = xyz[0];
        double y = xyz[1];
        double z = xyz[2];
        double Vx = VV[0];
        double Vy = VV[1];
        double Vz = VV[2];
        res[3] = ((-mu / Math.pow(r, 3)) * x - c_xa * (S_mid / (2 * m)) * ro * V * Vx +
                K * c_xa * (S_mid / (2 * m)) * ro * Math.pow(V, 2) * (m1 * Math.cos(gam) + l1 * Math.sin(gam))
                + x * Math.pow(omega, 2) + 2 * omega * Vy);

        res[4] = ((-mu / Math.pow(r, 3)) * y - c_xa * (S_mid / (2 * m)) * ro * V * Vy +
                K * c_xa * (S_mid / (2 * m)) * ro * Math.pow(V, 2) * (m2 * Math.cos(gam) + l2 * Math.sin(gam))
                + y * Math.pow(omega, 2) - 2 * omega * Vx);

        res[5] = ((-mu / Math.pow(r, 3)) * z - c_xa * (S_mid / (2 * m)) * ro * V * Vz +
                K * c_xa * (S_mid / (2 * m)) * ro * Math.pow(V, 2) * (m3 * Math.cos(gam) + l3 * Math.sin(gam)));
        return res;
    }
    public static double[] runge(double[] vect, double gam){
        double h=0.2;
        double[] res = new double[6];
        double[] k1 = right(vect, gam);
        double[] A = vekchi(k1, h * 0.5);
        double[] k2 = right(sum_ve(vect, A), gam);
        double[] k3 = right( sum_ve(vect, vekchi(k2, h * 0.5)), gam);
        double[] k4 = right(sum_ve(vect, vekchi(k3,h)), gam);

        double[] tmp = sum_ve(k1, k4);
        tmp = sum_ve(tmp, vekchi(k2,2));
        tmp = sum_ve(tmp, vekchi( k3,2));
        tmp = vekchi(tmp,h / 6 );
        res = sum_ve(vect, tmp);
        return res;
    }
    public static double[] VectProizv(double[] A, double[] B) {
        double[] Exit = new double[A.length];
        for (int i = 0; i < 3; i++) {
            Exit[i] = A[(i + 1) % 3] * B[(i + 2) % 3] - A[(i + 2) % 3] * B[(i + 1) % 3];
        }
        return Exit;
    }
    //Умножение вектора на число
    public static double[] vekchi(double A[], double B) {
        double[] C = new double[A.length];
        for (int i = 0; i < A.length; i++) {
            C[i] = A[i]*B;
        }
        return C;
    }
    //Сумма векторов
    public static double[] sum_ve(double A[], double B[]) {
        double[] C = new double[A.length];
        for (int i = 0; i < A.length; i++) {
            C[i] = A[i]+B[i];

        }
        return C;
    }
}