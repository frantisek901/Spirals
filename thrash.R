#### Thrash for ICA extended abstract






# Old graphs --------------------------------------------------------------



## Effects on size of components

ggplot(data, aes(x = size, fill = drawn, color = drawn)) +
  geom_density(alpha = 0.2) +
  labs(x = "Size of component as percentage of all agents", y = "Density",
       title = "Method of assigning boundary matters!",
       subtitle = "Random uniform distribution prefers components >80%\nassigning a constant value prefers components <10%.") +
  scale_x_continuous(breaks = seq(0, 100, 20), labels = function(x) paste0(x, "%")) +
  theme_minimal()
ggsave("Figs/Fig01-MethodOfAssigningBoundaryMatters.png", width = 6, height = 4)



data %>% mutate(neis = factor(neis)) %>%
  # filter(N_agents == 129) %>%
  ggplot(aes(x = size, fill = neis, color = neis)) +
  geom_density(alpha = 0.2) +
  labs(x = "Size of component as percentage of all agents", y = "Density",
       title = "Number of network neighbors matters!",
       subtitle = "The smaller neighborhood the more prefers smaller components <10%.") +
  scale_x_continuous(breaks = seq(0, 100, 20), labels = function(x) paste0(x, "%")) +
  theme_minimal()
ggsave("Figs/Fig02-SizeOfNeighborhoodMatters.png", width = 6, height = 4)



data %>% mutate(random_links = factor(random_links)) %>%
  # filter(neis == 4) %>%
  ggplot(aes(x = size, fill = random_links, color = random_links)) +
  geom_density(alpha = 0.2) +
  facet_grid(cols = vars(N_agents)) +
  labs(x = "Size of component as percentage of all agents", y = "Density",
       title = "Random re-wiring matters!",
       subtitle = "The smaller re-wiring more prefers smaller components <10%\n and number of agents in simulation intensifies this.") +
  scale_x_continuous(breaks = seq(20, 100, 40), labels = function(x) paste0(x, "%")) +
  theme_minimal()
ggsave("Figs/Fig03-RandomRewiringMatters.png", width = 6, height = 4)



##  Effects on length of simulation

data %>% filter(component_ID == 1) %>%  # Now we are interested in simulation as a whole so we have select just one compenent from whole simulation, since all components contain same information about whole simulation.
  ggplot(aes(x = step, fill = drawn, color = drawn)) +
  geom_density(alpha = 0.2) +
  labs(x = "Steps until the equilibrium", y = "Density",
       title = "Method of assigning boundary doesn't matter\nfor length of simulation.",
       #subtitle = "Random uniform distribution prefers components >80%\nassigning a constant value prefers components <10%."
  ) +
  scale_x_log10() + #breaks = seq(0, 100, 20), labels = function(x) paste0(x, "%")) +
  # scale_y_log10() +
  theme_minimal()
# NOTE: Since tehre is no effect we do not save the picture of graph.



data %>% filter(component_ID == 1 & step < 350) %>%  # Now we are interested in simulation as a whole so we have select just one compenent from whole simulation, since all components contain same information about whole simulation.
  mutate(neis = factor(neis)) %>%
  ggplot(aes(x = step, fill = neis, color = neis)) +
  geom_density(alpha = 0.2) +
  facet_grid(cols = vars(N_agents)) +
  labs(x = "Steps until the equilibrium", y = "Density",
       title = "Number of neighbors matters for length of simulation!",
       subtitle = "The more neighbors the shorter the simulation.\nIt's intensified by lower numbers of agents in simulation.",
       caption = "Note: for steps over 350 the neighborhood size doesn't matter,\nso we omit simulations that reach equilibrium after 350 steps."
  ) +
  # scale_x_log10() + #breaks = seq(0, 100, 20), labels = function(x) paste0(x, "%")) +
  # scale_y_log10() +
  theme_minimal()
ggsave("Figs/Fig04-SizeOfNeighborhoodMatters.png", width = 6, height = 4)



data %>% filter(component_ID == 1 & step < 350) %>%  # Now we are interested in simulation as a whole so we have select just one compenent from whole simulation, since all components contain same information about whole simulation.
  mutate(random_links = factor(random_links)) %>%
  ggplot(aes(x = step, fill = random_links, color = random_links)) +
  geom_density(alpha = 0.2) +
  facet_grid(cols = vars(N_agents)) +
  labs(x = "Steps until the equilibrium", y = "Density",
       title = "Random re-wiring matters for length of simulation!",
       subtitle = "The higher re-wiring the shorter the simulation.\nIt's intensified a bit by lower numbers of agents in simulation.",
       caption = "Note: for steps over 350 the neighborhood size doesn't matter,\nso we omit simulations that reach equilibrium after 350 steps."
  ) +
  # scale_x_log10() + #breaks = seq(0, 100, 20), labels = function(x) paste0(x, "%")) +
  # scale_y_log10() +
  theme_minimal()
ggsave("Figs/Fig05-RandomRewiringMatters.png", width = 6, height = 4)



## Effects on including simulation into main sample and reason of omitting

reason %>% mutate(neis = factor(neis)) %>%
  ggplot(aes(y = neis, fill = reason)) +
  geom_bar(position = "stack") +
  labs(title = "Size of neighborhood matters!",
       subtitle = "With larger neighborhood it's more likely\nthat simulation reaches equilibrium.\nBut constant is proportion of simulations\nwith all components smaller than 6.") +
  theme_minimal()
ggsave("Figs/Fig06-SizeOfNeighborhoodMatters.png", width = 6, height = 4)


reason %>% mutate(opinions = factor(opinions)) %>%
  ggplot(aes(y = opinions, fill = reason)) +
  geom_bar(position = "stack") +
  labs(title = "Number of opinions matters!",
       subtitle = "With more dimensions of opinion space is more likely\nthat all componets are of size smaller than 6.\nThe proportion of simulations not reaching equilibrium is more or less\nsame regardless the dimensions of opinion space.") +
  theme_minimal()
ggsave("Figs/Fig07-NumberOfOpinionsMatters.png", width = 6, height = 4)


reason %>% mutate(neis = factor(neis), opinions = factor(opinions)) %>%
  ggplot(aes(y = neis, fill = reason)) +
  facet_grid(vars(opinions)) +
  geom_bar(position = "fill") +
  scale_x_continuous(breaks = seq(0, 1, .20), labels = function(x) paste0(100*x, "%")) +
  labs(title = "Number of opinions and neighbors matter!\n('neis' on Y-axis, 'opinions' define panels)",
       subtitle = "Combination of larger neigborhood and less opinion dimensions supports\ninclusion of simulation into main file, i.e. simulation reaches equilibrium\nand at least some components are of size 6+. While more opinions\nfracture the public into groups smaller than 6, smaller neighborhood blocks\nreaching equilibrium.",
       x = "Proportion") +
  theme_minimal()
ggsave("Figs/Fig08-NumberOfOpinionsAndNeighborsMatter.png", width = 6, height = 4)




