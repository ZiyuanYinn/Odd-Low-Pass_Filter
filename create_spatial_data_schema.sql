CREATE TABLE spatial_data (
    id SERIAL PRIMARY KEY,
    l DOUBLE PRECISION NOT NULL,
    b DOUBLE PRECISION NOT NULL,
    parallax DOUBLE PRECISION,
    phot_g_mean_mag DOUBLE PRECISION,
    bp_rp DOUBLE PRECISION
);

CREATE INDEX idx_spatial_coordinates ON spatial_data (l, b);
CREATE INDEX idx_parallax ON spatial_data (parallax);


*.sql linguist-detectable=true
