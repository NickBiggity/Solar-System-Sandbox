/*///scr_simEngine()

    G = 6.67408 * power(10, 11);
    MS = 1.989 * power(10, 30);
    AU = 1.495978707 * power(10, 11);
    DAY = 86400;
    YEAR = DAY*365.25;
    
    bod = 11; //No. of bodies
    m = array_create(bod); //Array holding masses of bodies. Sun [0], Mercury [1]...
    inc = array_create(bod); //Array with inclinations of each body to the ecliptic
    X = array_create(bod); Y = array_create(bod); Z = array_create(bod); vx = array_create(bod); vy = array_create(bod); vz = array_create(bod); //Variables for X, Y positions and velocities in these directions.
    a[bod, 6] = 0; b[bod, 6] = 0; c[bod, 6] = 0; d[bod, 6] = 0; //RK variables for each body & variable we calculate, 0 = X, 1 = Y, 2 = Z, 3 = vx, 4 = vy, 5 = vz
    dt = 2500; tMax = YEAR*10; tYear = 0; t = 0; tTransit = 0;
    yprev = 0;
    
    //FILE * foutput;
    //foutput = fopen("output.txt", "w");
    
    //Masses of Sun and Planets
    m[0] = 1.989 * power(10, 30); m[1] = 3.285 * power(10, 23);// m[2] = 4.867E24; m[3] = 5.972e24; m[4] = 6.39E23;
    //m[5] = 1.8982E27; m[6] = 5.6834E26; m[7] = 8.68E25; m[8] = 1.0243E26; m[9] = 7.34767309E22;
    //m[10] = 0.95E21;
    
    //Inclinations of bodies
    inc[0] = 0; inc[1] = 7.005;// inc[2] = 3.395; inc[3] = 0; inc[4] = 1.850; inc[5] = 1.303;
    //inc[6] = 2.485; inc[7] = 0.773; inc[8] = 1.767; inc[9] = 0/*5.145; inc[10] = 10.593;
    
    X[0] = 0; Y[0] = 0; //Sun
    X[1] = -0.4*AU; Y[1] = 0; //Mercury
    /*X[2] = 0.7*AU; Y[2] = 0; //Venus
    X[3] = -AU; Y[3] = 0; //Earth
    X[4] = 1.5*AU; Y[4] = 0; //Mars
    X[5] = -5.2044*AU; Y[5] = 0; //Jupiter
    X[6] = 9.5826*AU; Y[6] = 0; //Saturn
    X[7] = -19.2184*AU; Y[7] = 0; //Uranus
    X[8] = 30.1*AU; Y[8] = 0; //Neptune
    X[9] = -AU*1.00256955529; Y[9] = 0; //Moon
    X[10] = 2.06*AU; Y[10] = 0; //Ceres
    
    for (i = 0; i < bod; i++) { //Sets correct Z & X positions for 3d (Also converts inclination values to radians)
            inc[i] = 0.0174532925*inc[i];
            X[i] = X[i]*cos(inc[i]);
            Z[i] = X[i]*sin(inc[i]);
    }
    
    vx[0] = 0; vy[0] = 0; //Sun at rest originally
    
    for (i = 1; i < bod; i++) { //Set initial velocities
        r = scr_pythag(X[i], Y[i], Z[i]);
        vx[i] = 0;
        vy[i] = sign(X[i])*sqrt((G*MS)/r);
        vz[i] = 0;
        vy[0] -= (m[i]*vy[i])/m[0]; //Conserves Momentum
    }
    vy[9] += sign(X[9])*sqrt((G*m[3])/point_distance_3d(X[3], Y[3], 0, X[9], Y[9], Z[9])); //Moon
    
    for (t = 0; t < tMax;) { //Time Loop
    for (n = 0; n < 100; t+=dt, n++) {
    yprev = 0;
    var sx, sy, sz; //Displacement between two bodies in X, Y and Z
    
    for (i = 0; i < bod; i++) { //Loops for a step in each body
    a[i][0] = vx[i];
    a[i][1] = vy[i];
    a[i][2] = vz[i];
    
                    for (j = 3; j < 6; j++) { //Resets velocity values for RK
                        a[i][j] = 0;
                        b[i][j] = 0;
                        c[i][j] = 0;
                        d[i][j] = 0;
                    }
    
    for (j = 0; j < bod; j++) {
    if (j != i) {
    sx = X[j]-X[i];
    sy = Y[j]-Y[i];
    sz = Z[j]-Z[i];
    var r = scr_pythag(sx, sy, sz);
    a[i][3] += ((G*m[j]*sx)/power(r, 3));
    a[i][4] += ((G*m[j]*sy)/power(r, 3));
    a[i][5] += ((G*m[j]*sz)/power(r, 3));
    }
    }
    }
    
    for (i = 0; i < bod; i++) {
    b[i][0] = vx[i] + (dt/2)*a[i][3];
    b[i][1] = vy[i] + (dt/2)*a[i][4];
    b[i][2] = vz[i] + (dt/2)*a[i][5];
    
    for (j = 0; j < bod; j++) {
    if (j != i) {
    sx = (X[j] + (dt/2)*a[j][0]) - (X[i] + (dt/2)*a[i][0]);
    sy = (Y[j] + (dt/2)*a[j][1]) - (Y[i] + (dt/2)*a[i][1]);
    sz = (Z[j] + (dt/2)*a[j][2]) - (Z[i] + (dt/2)*a[i][2]);
    var r = scr_pythag(sx, sy, sz);
    b[i][3] += ((G*m[j]*sx)/power(r, 3));
    b[i][4] += ((G*m[j]*sy)/power(r, 3));
    b[i][5] += ((G*m[j]*sz)/power(r, 3));
    }
    }
    }
    
    for (i = 0; i < bod; i++) {
    c[i][0] = vx[i] + (dt/2)*b[i][3];
    c[i][1] = vy[i] + (dt/2)*b[i][4];
    c[i][2] = vz[i] + (dt/2)*b[i][5];
    
    for (j = 0; j < bod; j++) {
    if (j != i) {
    sx = (X[j] + (dt/2)*b[j][0]) - (X[i] + (dt/2)*b[i][0]);
    sy = (Y[j] + (dt/2)*b[j][1]) - (Y[i] + (dt/2)*b[i][1]);
    sz = (Z[j] + (dt/2)*b[j][2]) - (Z[i] + (dt/2)*b[i][2]);
    var r = scr_pythag(sx, sy, sz);
    c[i][3] += ((G*m[j]*sx)/power(r, 3));
    c[i][4] += ((G*m[j]*sy)/power(r, 3));
    c[i][5] += ((G*m[j]*sz)/power(r, 3));
    }
    }
    }
    
    for (i = 0; i < bod; i++) {
    d[i][0] = vx[i] + dt*c[i][3];
    d[i][1] = vy[i] + dt*c[i][4];
    d[i][2] = vz[i] + dt*c[i][5];
    
    for (j = 0; j < bod; j++) {
    if (j != i) {
    sx = (X[j] + dt*c[j][0]) - (X[i] + dt*c[i][0]);
    sy = (Y[j] + dt*c[j][1]) - (Y[i] + dt*c[i][1]);
    sz = (Z[j] + dt*c[j][2]) - (Z[i] + dt*c[i][2]);
    var r = scr_pythag(sx, sy, sz);
    d[i][3] += ((G*m[j]*sx)/power(r, 3));
    d[i][4] += ((G*m[j]*sy)/power(r, 3));
    d[i][5] += ((G*m[j]*sz)/power(r, 3));
    }
    }
    }
    for (i = 0; i < bod; i++) {
    X[i] += (dt/6)*(a[i][0] + 2*b[i][0] + 2*c[i][0] + d[i][0]);
    Y[i] += (dt/6)*(a[i][1] + 2*b[i][1] + 2*c[i][1] + d[i][1]);
    Z[i] += (dt/6)*(a[i][2] + 2*b[i][2] + 2*c[i][2] + d[i][2]);
    vx[i] += (dt/6)*(a[i][3] + 2*b[i][3] + 2*c[i][3] + d[i][3]);
    vy[i] += (dt/6)*(a[i][4] + 2*b[i][4] + 2*c[i][4] + d[i][4]);
    vz[i] += (dt/6)*(a[i][5] + 2*b[i][5] + 2*c[i][5] + d[i][5]);
    }
    
    if ((Y[3]*yprev < 0) && (yprev > 0)) {
    tTransit = ((t - dt) - 0) + yprev/vy[3] - tYear;
    //printf("%lf\n", tTransit/DAY);
    tYear += tTransit;
    }
    }
    for (i = 0; i < bod; i++) {
    //fprintf(foutput, "%lf\t%lf\t%lf\t", X[i]/AU, Z[i]/AU, Y[i]/AU);
    }
            //fprintf(foutput, "%lf\n", 0.5*m[10]*pythag(vx[10], vy[10])*pythag(vx[10], vy[10]));
    //fprintf(foutput, "\n");
            //printf("%lf%%\n", (t/tMax)*100);
    
    }

//fclose(foutput);*/
