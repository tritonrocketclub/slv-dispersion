/**
*This method assumes that the piece of debris is always falling at terminal velocity.
*The time required to fall in an altitude band defined by hi to hi+1 is called the dwell time.
*It computes the dwell time by dividing the differential altitude (hi - hi+1) by the 
*terminal velocity, VT in that altitude range, where VT = sqrt(2B/rho).
*It computes VT by integrating over the altitudeinterval with delta-t
* changing with altitude. It then assumes that the fragment moves
*horizontally exactly at the speed of the wind in the altitude interval during the time that it
*is falling through the altitude interval.
*/

public double[] Drift(double B, double h, double[] ewind, double[] nwind)
    {
        
        final double k1 = 0.0065/288.15;				//   L/To
        final double k2 = (9.8*0.029)/(8.31*0.0065);			//   (gM)/(RL)
        final double k3 = 10*Math.sqrt((5*0.029*101.325)/(B*8.31));	//   10sqrt((5MPo)/BR)
        final int INCREMENT = 1000;					//   increment in elevation (hi - hi-1)
        final int LOOPS = 10;						//   iterations of RRAM during integration
        
        double dE = 0;							//   east-west displacement due to wind
        double dN = 0;							//   north-south displacement due to wind
        int i = (int)(h+INCREMENT/2)/INCREMENT-1;			//   rounds to nearest multiple of the increment
		
        while(i>=0)
        {
            double t = 0;						//   dwell time spent in this altitude band of atmosphere
            double dh = INCREMENT/LOOPS;				//   height step in integration through the band
            for(int n = 1; n <= LOOPS; n++)				//   RRAM integration to calculate dwell time
            {
		double height = (i-1)*(INCREMENT)+n*dh;
                t = t + Math.sqrt((Math.pow(1-k1*height, k2))/(288.15-0.0065*height))*dh;
            }
            t = Math.abs(k3*t);
            dE = dE + t*ewind[i];		 
            dN = dN + t*nwind[i];
            i--;
        }
        double[] result = {dE,dN};
        return result;
    }
