
function RASTER = smoothRaster ( RASTER , convKernel )

if iscell(RASTER),  %if RASTER is a cell returned by getRaster
    
    for g = 1 : length(RASTER),
        for r = 1 : size(RASTER{g},1),
                %since getRaster allows exclusion of parts of the raster, we
                %need to find valid chuncks within the raster before the
                %convolution, otherwise the nans in the raster will corrupt the
                %result of the convolution
            I = find(~isnan(RASTER{g}(r,:)));
            start_ind = I(diff([-Inf I])~=1);
            end_ind = I(diff([I Inf])~=1);
                %do the convolution
            for s = 1 : length(start_ind),
                RASTER{g}(r,start_ind(s):end_ind(s)) = nanconv(RASTER{g}(r,start_ind(s):end_ind(s)),convKernel,'same');
            end;
        end;
    end;
    
else                %if RASTER is a matrix with rows corresponding to trials
    
        for r = 1 : size(RASTER,1),
                %since getRaster allows exclusion of parts of the raster, we
                %need to find valid chuncks within the raster before the
                %convolution, otherwise the nans in the raster will corrupt the
                %result of the convolution
            I = find(~isnan(RASTER(r,:)));
            start_ind = I(diff([-Inf I])~=1);
            end_ind = I(diff([I Inf])~=1);
                %do the convolution
            for s = 1 : length(start_ind),
                RASTER(r,start_ind(s):end_ind(s)) = nanconv(RASTER(r,start_ind(s):end_ind(s)),convKernel,'same');
            end;
        end;
    
end;





