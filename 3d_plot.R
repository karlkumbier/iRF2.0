load('intermediate_results.Rdata')
# plot 3 features (wt_ZLD, gt2, kr1) and the color is the label. 
#   It can be seen that when gt2 or ZLD or kr is small, the label is dominantly zero.
library(plotly)
data = as.data.frame(cbind(X,Y))
colnames(data)
p3 <- plot_ly(data[train.id[indices],], x = ~wt_ZLD, y = ~gt2, z = ~kr1,
             marker = list(color = ~Y, opacity = .5, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'ZLD'),
                      yaxis = list(title = 'gt'),
                      zaxis = list(title = 'kr')))
p3

# ZLD v.s. gt
p2a <- plot_ly(data[train.id, ], x = ~wt_ZLD, y = ~gt2, 
              marker = list(color = ~Y, opacity = .5, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'ZLD'),
                      yaxis = list(title = 'gt')))
p2a


# gt v.s. kr
p2b <- plot_ly(data[train.id[indices], ], x = ~gt2, y = ~kr1,
               marker = list(color = ~Y, opacity = .5, colorscale = 'RdBu', showscale = TRUE)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'gt2'),
                      yaxis = list(title = 'kr')))
p2b

# gt v.s. D1
p2c <- plot_ly(data[train.id, ], x = ~gt2, y = ~D1,
               marker = list(color = ~Y, opacity = .5, colorscale='RdBu', showscale = TRUE)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'gt2'),
                      yaxis = list(title = 'D1')))
p2c

data %>% filter(1:nrow(data) %in% train.id) %>% select(gt2, D1, Y) %>% 
  mutate(gt2 = gt2 > 1) %>% mutate(D1 = D1 > .5) %>% 
  group_by(D1, gt2) %>% summarise(mean(Y))
