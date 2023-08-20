function show_distance(lat1, lon1, lat2, lon2)
    % Coordinates of two locations (latitude, longitude)    
    distance = haversine_distance(lat1, lon1, lat2, lon2);
    fprintf('Distance: %.2f km\n', distance);
end