function hex = rgb2hex(rgb)
    % Convert rgb color to hexdecimal

    if max(rgb(:))<=1
        rgb = round(rgb*255);
    else
        rgb = round(rgb);
    end

    hex(:,2:7) = reshape(sprintf('%02X',rgb.'),6,[]).';
    hex(:,1) = '#';
end