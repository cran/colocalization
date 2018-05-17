########
## 2D ##
########
set.seed(1234)
x <- runif(300, min = -1, max = 1)
y <- runif(300, min = -1, max = 1)
red <- data.frame(x,y, color = "red")
x <- runif(50, min = -1, max = 1)
y <- runif(50, min = -1, max = 1)
green <- data.frame(x,y, color = "green")
mydata <- rbind(red,green)

# Class of the returned results
#result_z <- nsinc.z(data = mydata, membership = "color", dim = 2)
#expect_that(result_z, is_a("colocal"))

# Warnings
#expect_that(nsinc.z(data = mydata, membership = "color", dim = 2,
#                             r.model = "full", r.adjust = 0),
#            gives_warning("The r.adjust must be a positive number no larger than half of the difference between r.max and r.min!"))
#expect_that(nsinc.z(data = mydata, membership = "color", dim = 2,
#                             r.model = "full", r.adjust = -1),
#            gives_warning("The r.adjust must be a positive number no larger than half of the difference between r.max and r.min!"))

distance.matrix <- as.matrix(dist(mydata[,c("x","y")]))
r.min.post.data <-  min(distance.matrix[distance.matrix!=0])
r.max.post.data <- max(distance.matrix[distance.matrix!=0])

#expect_that(nsinc.z(data = mydata, membership = "color", dim = 2,
#                             r.model = "full", r.adjust = ((r.max.post.data-r.min.post.data)/2 + 0.5)),
#            gives_warning("The r.adjust must be a positive number no larger than half of the difference between r.max and r.min!"))

#expect_that(nsinc.z(data = mydata, membership = "color", dim = 2,
#                             r.model = "full", r.count = 0),
#            gives_warning("r.count must be a positive integer"))
#expect_that(nsinc.z(data = mydata, membership = "color", dim = 2,
#                             r.model = "full", r.count = -1),
#            gives_warning("r.count must be a positive integer"))

#expect_that(nsinc.z(data = mydata, membership = "color", dim = 2,
#                             r.model = "other", r.min = 0.01, r.max = 0.5, r.count = 0),
#            gives_warning("r.count must be a positive integer"))
#expect_that(nsinc.z(data = mydata, membership = "color", dim = 2,
#                             r.model = "other", r.min = 0.01, r.max = 0.5, r.count = -1),
#            gives_warning("r.count must be a positive integer"))


# Errors
expect_that(nsinc.z(data = mydata, membership = "member", dim = 2),
            throws_error("There is no column names in the data called member!"))

expect_that(nsinc.z(data = mydata, membership = "color", dim = 4),
            throws_error("dim must be either 2 or 3!"))

colnames(mydata) <- c("a","y","color")
expect_that(nsinc.z(data = mydata, membership = "color", dim = 2),
            throws_error("Data must contain a 'x' column!"))

colnames(mydata) <- c("x","b","color")
expect_that(nsinc.z(data = mydata, membership = "color", dim = 2),
            throws_error("Data must contain a 'y' column!"))

colnames(mydata) <- c("x", "y", "color")

box <- data.frame(min=-1, xmax=1, ymin=-1, ymax=1)
expect_that(nsinc.z(data = mydata, membership = "color", dim = 2, box = box),
            throws_error("'box' must be a dataframe containing columns 'xmin','xmax','ymin' and 'ymax'!"))

box <- data.frame(xmin=2, xmax=1, ymin=-1, ymax=1)
expect_that(nsinc.z(data = mydata, membership = "color", dim = 2, box = box),
            throws_error("'xmax' or 'ymax' must be larger than 'xmin' or 'ymin' in 'box'!"))

expect_that(nsinc.z(data = red, membership = "color", dim = 2),
            throws_error("There must be at least two memberships of signals in the input data!"))

box <- data.frame(xmin=2, xmax=3, ymin=2, ymax=3)
expect_that(nsinc.z(data = mydata, membership = "color", dim = 2, box = box),
            throws_error("There must be at least two memberships of signals enclosed in the study region!"))

expect_that(nsinc.z(data = mydata, membership = "color", dim = 2,
                             r.model = "other", r.max = 0.5),
            throws_error("If choose the 'other' for r.model, then r.min must be specified by the user!"))
expect_that(nsinc.z(data = mydata, membership = "color", dim = 2,
                             r.model = "other", r.min = 0.01),
            throws_error("If choose the 'other' for r.model, then r.max must be specified by the user!"))

expect_that(nsinc.z(data = mydata, membership = "color", dim = 2,
                             r.model = "other", r.min = r.min.post.data-0.1, r.max = 0.5),
            throws_error("r.min must be between the smallest and half of the largest interpoint distances"))
expect_that(nsinc.z(data = mydata, membership = "color", dim = 2,
                             r.model = "other", r.min = 0.5*r.max.post.data+0.1, r.max = 0.5),
            throws_error("r.min must be between the smallest and half of the largest interpoint distances"))

expect_that(nsinc.z(data = mydata, membership = "color", dim = 2,
                             r.model = "other", r.min = 0.01, r.max = r.min.post.data-0.1),
            throws_error("r.max must be between the smallest and half of the largest interpoint distances"))
expect_that(nsinc.z(data = mydata, membership = "color", dim = 2,
                             r.model = "other", r.min = 0.01, r.max = 0.5*r.max.post.data+0.1),
            throws_error("r.max must be between the smallest and half of the largest interpoint distances"))

expect_that(nsinc.z(data = mydata, membership = "color", dim = 2,
                             r.model = "other", r.min = 0.5, r.max = 0.01),
            throws_error("The r.min must be smaller than r.max!"))

expect_that(nsinc.z(data = mydata, membership = "color", dim = 2,
                    r.model = "other", r.min = 0.01, r.max = 0.5, r.adjust = -1),
            throws_error("The r.adjust must be a nonnegative number smaller than half of the difference between r.max and r.min!"))
expect_that(nsinc.z(data = mydata, membership = "color", dim = 2,
                    r.model = "other", r.min = 0.01, r.max = 0.5, r.adjust = 0.25),
            throws_error("The r.adjust must be a nonnegative number smaller than half of the difference between r.max and r.min!"))

expect_that(nsinc.z(data = mydata, membership = "color", dim = 2, r.model = "Bayesian"),
            throws_error("r.model must be one of 'full'"))

expect_that(nsinc.z(data = mydata, membership = "color", dim = 2,
                             strata = TRUE, base.member = "blue"),
            throws_error("The specified base membership 'blue' is not found in the provided data!"))
