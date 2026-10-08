test_that("moss_data runs without error", {
      expect_no_error(moss())
      expect_no_error(moss("poly"))
})

test_that("moss companion layers load in raster and polygon formats", {
      ps <- moss()
      for(d in c("protection", "popdens")){
            r <- moss(data = d)
            expect_s4_class(r, "SpatRaster")
            expect_equal(names(r), d)
            expect_true(terra::compareGeom(r, ps$spatial, stopOnError = FALSE))

            p <- moss("polygon", data = d)
            expect_s3_class(p, "sf")
            expect_true(d %in% names(p))
            expect_equal(nrow(p), ps$n_sites)
            expect_equal(p[[d]], as.vector(terra::values(r)))
      }
      expect_error(moss(data = "elevation"))
})

test_that("moss companion layers are valid prioritization inputs", {
      ps <- moss()
      v <- site_values(moss(data = "protection"))[ps$occupied]
      expect_true(all(is.finite(v)) && min(v) >= 0 && max(v) <= 1)
      v <- site_values(moss(data = "popdens"))[ps$occupied]
      expect_true(all(is.finite(v)) && min(v) >= 0)
})
