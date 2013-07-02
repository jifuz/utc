function rgb=visualize_combinePriorIm(prior,sourceImage)

%0->unknown
% label('1')='parking';
% label('2')='footway';
% label('3')='highway';
% label('4')='building';
% label('5')='transportation_objects';

% 0 sky 128 128 128
% 1 tree 128 128 0
% 2 asphalt_floor 128 64 128
% 3 grass 0 128 0
% 4 building 128 0 0
% 5 object 192 128 0
% 6 concrete_floor 0 128 128
% 7 gravel_floor 199 21 133
% 8 vehicle 0 0 128
labelData=prior;

  [m,n]=size(labelData);
    red = [128 199 0 128 128 192 ];%  128 128 0 128 192 0 199 0];
    green = [128 21 128 64 0  128 ];%128 64 128 0 128 128 21 0];
    blue = [128 133 128 128 0 0 0 ];%128 0 0 0 128 133 128];
    im=double(sourceImage);
    rgb = zeros(m,n,3);
    for i = 1:m
        for j = 1:n
            if (labelData(i,j) >= 0)
                rgb(i,j,1) = 0.5*(red(labelData(i,j)+1)+im(i,j,1))/255.0;
                rgb(i,j,2) = 0.5*(green(labelData(i,j)+1)+im(i,j,2))/255.0;
                rgb(i,j,3) = 0.5*(blue(labelData(i,j)+1)+im(i,j,3))/255.0;
            end
        end
    end
%     figure;
%     imagesc(rgb);

end