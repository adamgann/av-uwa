%AV Auxiliary-vector sequence generation for adaptive interference
%   suppressive filtering.
%   [g, w]=AV(R, v,c, t, maxiters) produces the matrix g, the columns of
%   which are the length-D auxiliary vectors g(:,1), g(:,2), ..., g(:,N),
%   and the matrix w, the columns of which are the sequence of size N-by-1 
%   av filters w(:,1), w(:,2), ..., w(:,N).
%   The converging point of the sequence w(:,1), w(:,2), ..., w(:,N)is 
%   winfty=c*inv(R)*v/(v'*inv(R)*v).
%   For each filter of the sequence, say w(:,n), it holds w(:,n)'*v=c.
%
%   If t=1, the sequence terminates upon convergence; i.e., 
%   when w(:,N)=w(:,N-1).
%   If t=0, the sequence terminates at w(:,N) for N=maxiters,
%   possibly prior to convergence.
%
%   A note on filter design: 
%   If R is the true received signal autocorrelation matrix and v is the 
%   input/desired-output cross-correlation vector (i.e., the 
%   mathced filter), the converging point of the sequence is the 
%   MMSE/max-SINR filter. 
%   If R and v are estimates upon M observations, the converging point
%   of the sequence is the corresponding SMI. In this case, non-asymptotic
%   elements of the AV filter sequence exhibit higher SINR and lower
%   MVDR-filter-estimation error than the converging point, when M is
%   relatively (w.r.t. D) small. 
%
%   Author:
%   Panos P. Markopoulos
%   Ph.D. Student and Research Assistant
%   Signals, Communications, and Networking Research Group
%   Electrical Engineering Department
%   University at Buffalo, The State University of New York
%   Email: pmarkopo@buffalo.edu or pados@buffalo.edu
%   Web: http://www.acsu.buffalo.edu/~pmarkopo/
%
%   Some references: 
%   1) D. A. Pados and S. N. Batalama, ``Low-complexity blind detection
%   of DS/CDMA signals: Auxiliary-vector receivers," IEEE Transactions on
%   Communications, vol. 45, pp. 1586-1594, Dec. 1997.
%   2) D. A. Pados, F. J. Lombardo, and S. N. Batalama,``Auxiliary-vector
%   filters and adaptive steering for DS/CDMA single-user detection," IEEE
%   Transactions on Vehicular Technology, vol. 48, pp. 1831-1839, 
%   Nov. 1999.
%   3)D. A. Pados and S. N. Batalama, ``Joint space-time auxiliary-vector
%   filtering for DS/CDMA systems with antenna arrays," IEEE Transactions
%   on Communications, vol. 47, pp. 1406- 1415, Sept. 1999.
%   4) D. A. Pados and G. N. Karystinos, ``An iterative algorithm for the
%   computation of the MVDR filter," IEEE Transactions on Signal
%   Processing, vol. 49, pp. 290-300, Feb. 2001.
%   5) R. Grover, D. A. Pados, and M. J. Medley, ``Subspace direction
%   finding with an Auxiliary-Vector basis," IEEE Transactions on Signal
%   Processing, vol. 55, pp. 758-763, Feb. 2007.
%   6) P. P. Markopoulos, S. Kundu, and D. A. Pados, ``Small-sample-support
%   suppression of interference to PN-masked data," IEEE Transactions on
%   Communications, vol. 61, pp. 2979-2987, July 2013.
%   7) P. P. Markopoulos, S. Kundu, and D. A. Pados, �Short-data-record
%   filtering of PN-masked data," in Proceedings 38th IEEE International
%   Conference on Acoustics, Speech, and Signal Processing (ICASSP 2013),
%   Vancouver, Canada, May 2013, pp. 4559-4563.
%
%   **Inquiries regarding the script provided below are cordially welcome. 
%   In case you spot a bug, please let me know. 
%   If you use some piece of code for your own work, please cite 
%   appropriately the articles above.** 

 
 
function [g, w]=av(R, v,c, t, maxiters)

w0=c*v/norm(v)^2;
epsilon=0.0005;

P=(eye(length(w0))-(1/norm(w0)^2)*(w0*w0'));
w(:,1)=w0;
for n=2:maxiters     
    gn=P*R*w(:,n-1);    
    if (norm(gn)<=epsilon) && t
        break;
    else
        g(:,n-1)=gn;       
        mn=(1/(gn'*R*gn))*(gn'*R*w(:,n-1));
        w(:,n)=w(:,n-1)-mn*gn;        
    end
end

