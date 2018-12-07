function [ Dice ] = DiceSimilarity3D( nx, ny,nz, II0, II1 ,LabelSize)
%DiceSimilarity Compute DiceSimilarity between given moving image II0 and fixed
%image II1 . 
%Input: nx,ny, nz = size of image in x, y, and z direction 
%        II0 = moving image
%        II1 = fixed image
%        LabelSize = number of label 
%Output  Dice =  Dice similarity of each label between image II0 and II1        

I0=reshape(II0, [ny nx nz]);
I1=reshape(II1, [ny nx nz]);


Dice=zeros(1,LabelSize);

for Label=0:LabelSize-1

  C1=numel(find(I0==Label));
            C2=numel(find(I1==Label));
            count=0;
            for i=1:nx
                for j=1:ny
                    for k =1:nz
                        if I0(j,i,k)==(Label) && I1(j,i,k)==(Label)
                            count=count+1;
                        end 
                    end
                end   
            end 
            Dice(Label+1)=(2*count)/(C1+C2);
end

