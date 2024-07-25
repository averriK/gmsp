function [PSA, PSV, SD, SA, SV, OUT] = responseSpectra(xi, sPeriod, gacc, dt)
%   PSEUDO SPECTRAL ACCELERATION, PSEUDO SPECTRAL VELOCITY AND SPECTRAL
%   DISPLACEMENT COMPUTATION
%
%   Syntax:
%     [PSA, PSV, SD, SA, SV, OUT] = responseSpectra(xi, sPeriod, gacc, dt)
%
%   Input:
%           xi = ratio of critical damping (e.g., 0.05)
%      sPeriod = vector of spectral periods
%         gacc = input acceleration time series in cm/s�
%           dt = sampling interval in seconds (e.g., 0.005)
%
%   Output:
%          PSA = Pseudo-spectral acceleration ordinates
%
%          PSV = Pseudo-spectral velocity ordinates
%
%           SD = Spectral displacement ordinates
%
%           SA = Spectral acceleration ordinates
%
%           SV = Spectral velocity ordinates
%
%          OUT = Time series of acceleration, velocity and displacemet response of SDF
%          oscillator
%
%   Example: Use demo.m in the zip file or copy paste commands below
%
%             load('gacc.mat');
%             dt = 0.005;
%             xi = 0.05;
%        sPeriod = [0.01,0.02,0.022,0.025,0.029,0.03,0.032,0.035,0.036,...
%                   0.04,0.042,0.044,0.045,0.046,0.048,0.05,0.055,0.06,0.065,0.067,0.07,...
%                   0.075,0.08,0.085,0.09,0.095,0.1,0.11,0.12,0.125,0.13,0.133,0.14,0.15,...
%                   0.16,0.17,0.18,0.19,0.2,0.22,0.24,0.25,0.26,0.28,0.29,0.3,0.32,0.34,...
%                   0.35,0.36,0.38,0.4,0.42,0.44,0.45,0.46,0.48,0.5,0.55,0.6,0.65,0.667,...
%                   0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,...
%                   2,2.2,2.4,2.5,2.6,2.8,3,3.2,3.4,3.5,3.6,3.8,4,4.2,4.4,4.6,4.8,5,7.5,10];
%
%   [PSA, PSV, SD] = responseSpectra(xi, sPeriod, gacc, dt);
%
%   Reference:
%
%   Wang, L.J. (1996). Processing of near-field earthquake accelerograms:
%   Pasadena, California Institute of Technology.
%
%   Acknowledgement: 
%
%   In plotSpectra function, I benefitted from an improved subplot
%   function (subplot1) by Eran Ofek
% 
%   THIS SOFTWARE IS PROVIDED "AS IS" AND ANY
%   EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
%   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER BE LIABLE
%   FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
%   BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
%   WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
%   OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
%   ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%   Written by Dr. Erol Kalkan, P.E. (ekalkan@usgs.gov)
%   $Revision: 3.0 $  $Date: 2017/03/27 17:03:00 $

vel = cumtrapz(gacc)*dt;
disp = cumtrapz(vel)*dt;

% Spectral solution
for i = 1:length(sPeriod)
    omegan = 2*pi/sPeriod(i);
    C = 2*xi*omegan;
    K = omegan^2;
    y(:,1) = [0;0];
    A = [0 1; -K -C]; Ae = expm(A*dt); AeB = A\(Ae-eye(2))*[0;1];
    
    for k = 2:numel(gacc)
        y(:,k) = Ae*y(:,k-1) + AeB*gacc(k);
    end
    
    displ = (y(1,:))';                          % Relative displacement vector (cm)
    veloc = (y(2,:))';                          % Relative velocity (cm/s)
    foverm = omegan^2*displ;                    % Lateral resisting force over mass (cm/s�)
    absacc = -2*xi*omegan*veloc-foverm;         % Absolute acceleration from equilibrium (cm/s�)
    
    % Extract peak values
    displ_max(i) = max(abs(displ));             % Spectral relative displacement (cm)
    veloc_max(i) = max(abs(veloc));             % Spectral relative velocity (cm/s)
    absacc_max(i) = max(abs(absacc));           % Spectral absolute acceleration (cm/s�)
    
    foverm_max(i) = max(abs(foverm));           % Spectral value of lateral resisting force over mass (cm/s�)
    pseudo_acc_max(i) = displ_max(i)*omegan^2;  % Pseudo spectral acceleration (cm/s)
    pseudo_veloc_max(i) = displ_max(i)*omegan;  % Pseudo spectral velocity (cm/s)
    
    PSA(i) = pseudo_acc_max(i);                 % PSA (cm/s�)
    SA(i)  = absacc_max(i);                     % SA (cm/s�)
    PSV(i) = pseudo_veloc_max(i);               % PSV (cm/s)
    SV(i)  = veloc_max(i);                      % SV (cm/s)
    SD(i)  = displ_max(i);                      % SD  (cm)
   
    % Time series of acceleration, velocity and displacement response of
    % SDF oscillator
    OUT.acc(:,i) = absacc;
    OUT.vel(:,i) = veloc;
    OUT.disp(:,i) = displ;
end