function [PHI,PHIU,PHIV,PHIW] = GaussPhi_fb(ac,BIGX, BIGY, BIGZ, Jm,knotvectorU,knotvectorV,knotvectorW, Coeff,param)
%#codegen
% This function computes the basis function (phi) along with the first
% derivative in each parametric direction (phi_u, phi_v, phi_w) at each
% gauss point in the active element

% This function also stores the coordinates of the gauss points for the intensities interpolated at the
% gauss points in each active cell stored in BIGX, BIGY, BIGZ

% H matrix stores the (size of each active element/2) at each active element

% INPUT:
% ac: active element array
% Em: element array struct form
% knotvectorU: knot vector in u direction
% knotvectorV: knot vector in v direction
% knotvectorW: knot vector in w direction
% Coeff: coefficient matrix for each active element
% param: the struct variable storing the parameters

% OUTPUT:
% PHI: cell array containing the basis functions values at each gauss point in each active
% element
% PHIU: cell array containing the derivative in u direction
%basis functions values at each gauss point in each active element
% PHIV: cell array containing the derivative in v direction basis functions values at each gauss point in each active
% element
% PHIW: cell array containing the derivative in w direction basis functions values at each gauss point in each active
% element
% BIGX, BIGY, BIGZ: the coordinates of the physical location of the gauss
% points in each active element
% H: the size of each active element divided by 2 (Jacobian of the transformation)
%%
% order of gaussian quadrature
orderGauss = param.orderGauss;

%degree of splines
pU = param.pU;
pV = param.pV;
pW = param.pW;

%number of splines in each direction
nobU = param.nobU;
nobV = param.nobV;
nobW = param.nobW;
ac_ct = size(ac,1);

%the gauss points of the corresponding gaussian quadrature order
%[si1,mmm] = ggquad(orderGauss);
%[si2,mmm] = ggquad(orderGauss);
%[si3,mmm] = ggquad(orderGauss);

PHI =  cell(ac_ct,1);
PHIU = cell(ac_ct,1);
PHIV = cell(ac_ct,1);
PHIW = cell(ac_ct,1);

for i = 1:ac_ct
    
    %cell_index = ac(i,1);
    cell_level = ac(i,2);
    
    SBX = BIGX(1+(i-1)*orderGauss:i*orderGauss,1:orderGauss,1:orderGauss);
    SBY = BIGY(1+(i-1)*orderGauss:i*orderGauss,1:orderGauss,1:orderGauss);
    SBZ = BIGZ(1+(i-1)*orderGauss:i*orderGauss,1:orderGauss,1:orderGauss);
    
    %coefficient matrix of non zero splines
    SB = Jm(i).nzsplines;
    supp_size = size(SB,1);
    gg_coeff = Coeff(i).mat;
    
    %knot vectors at the refinement level of the element
    knotU = knotvectorU{cell_level,1};
    knotV = knotvectorV{cell_level,1};
    knotW = knotvectorW{cell_level,1};
    unobg = nobU(cell_level,1);
    vnobg = nobV(cell_level,1);
    wnobg = nobW(cell_level,1);
    
    g_phi = zeros(supp_size,orderGauss,orderGauss);
    g_phiu = zeros(supp_size,orderGauss,orderGauss);
    g_phiv = zeros(supp_size,orderGauss,orderGauss);
    g_phiw = zeros(supp_size,orderGauss,orderGauss);
    
    %store the phis in the corresponding cell array
    for gg1 = 1:orderGauss
        for gg2 = 1:orderGauss
            for gg3 = 1:orderGauss
                uu = SBX(gg3,gg2,gg1);
                vv = SBY(gg3,gg2,gg1);
                ww = SBZ(gg3,gg2,gg1);
                
                uu = max(0,uu);
                uu = min(knotU(1,end),uu);
                vv = max(0,vv);
                vv = min(knotV(1,end),vv);
                ww = max(0,ww);
                ww = min(knotW(1,end),ww);
                
                uknot = FindSpan(unobg-1,pU,uu,knotU) + 1;
                vknot = FindSpan(vnobg-1,pV,vv,knotV) + 1;
                wknot = FindSpan(wnobg-1,pW,ww,knotW) + 1;

                RRD1 = Der1BasisFun(uknot-1,uu,pU,knotU);
                RRD2 = Der1BasisFun(vknot-1,vv,pV,knotV);
                RRD3 = Der1BasisFun(wknot-1,ww,pW,knotW);
                
                RRD = zeros((pU+1),(pV+1),(pW+1));
                RRDU = zeros((pU+1),(pV+1),(pW+1));
                RRDV = zeros((pU+1),(pV+1),(pW+1));
                RRDW = zeros((pU+1),(pV+1),(pW+1));
                
                for m3 = 1:(pU+1),
                    for m2 = 1:(pV+1),
                        for m1 = 1:(pW+1),
                            RRD(m1,m2,m3) = RRD1(1,m1)*RRD2(1,m2)*RRD3(1,m3);
                            RRDU(m1,m2,m3) = RRD1(2,m1)*RRD2(1,m2)*RRD3(1,m3);
                            RRDV(m1,m2,m3) = RRD1(1,m1)*RRD2(2,m2)*RRD3(1,m3);
                            RRDW(m1,m2,m3) = RRD1(1,m1)*RRD2(1,m2)*RRD3(2,m3);
                        end
                    end
                end
                
                inc=0;
                phii = zeros((pU+1)*(pV+1)*(pW+1),1);
                phiiu = zeros((pU+1)*(pV+1)*(pW+1),1);
                phiiv = zeros((pU+1)*(pV+1)*(pW+1),1);
                phiiw = zeros((pU+1)*(pV+1)*(pW+1),1);
                
                for m1 = 1:size(RRD,1),
                    for m2=1:size(RRD,2),
                        for m3 = 1:size(RRD,3),
                            phii((pU+1)*(pV+1)*(pW+1)-inc,1)= RRD(m3,m2,m1);
                            phiiu((pU+1)*(pV+1)*(pW+1)-inc,1)= RRDU(m3,m2,m1);
                            phiiv((pU+1)*(pV+1)*(pW+1)-inc,1)= RRDV(m3,m2,m1);
                            phiiw((pU+1)*(pV+1)*(pW+1)-inc,1)= RRDW(m3,m2,m1);
                            inc=inc+1;
                        end
                    end
                end
                phi_pi = gg_coeff*phii;
                phi_piu =gg_coeff*phiiu;
                phi_piv = gg_coeff*phiiv;
                phi_piw = gg_coeff*phiiw;
                g_phi(:,gg3,gg2,gg1) = phi_pi;
                g_phiu(:,gg3,gg2,gg1) = phi_piu;
                g_phiv(:,gg3,gg2,gg1) = phi_piv;
                g_phiw(:,gg3,gg2,gg1) = phi_piw;
            end
        end
    end
    
    
    PHI{i} = single(g_phi);
    PHIU{i} = single(g_phiu);
    PHIV{i} = single(g_phiv);
    PHIW{i} = single(g_phiw);
    
end



end