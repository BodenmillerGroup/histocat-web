<template>
  <v-expansion-panels :value="0">
    <v-expansion-panel>
      <v-expansion-panel-header disable-icon-rotate>
        <template v-slot:actions>
          <v-tooltip bottom>
            <template v-slot:activator="{ on }">
              <v-btn icon v-on="on" @click="deleteStep(step)">
                <v-icon color="secondary">mdi-delete-outline</v-icon>
              </v-btn>
            </template>
            <span>Delete step</span>
          </v-tooltip>
        </template>
        Neighborhood graph
      </v-expansion-panel-header>
      <v-expansion-panel-content>
        <v-card flat>
          <v-row dense>
            <v-col>
              <v-text-field
                label="Number of neighbors"
                hint="The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation."
                v-model.number="step.nNeighbors"
                type="number"
                min="2"
                max="100"
                step="1"
                :rules="[required]"
              />
            </v-col>
            <v-col>
              <v-select label="Metric" :items="metrics" v-model="step.metric" />
            </v-col>
            <v-col>
              <v-text-field label="Random seed" v-model.number="step.randomState" type="number" min="0" step="1" />
            </v-col>
          </v-row>
        </v-card>
      </v-expansion-panel-content>
    </v-expansion-panel>
  </v-expansion-panels>
</template>

<script lang="ts">
import { Component, Prop, Vue } from "vue-property-decorator";
import { required } from "@/utils/validators";

@Component
export default class NeighborsStepEditor extends Vue {
  @Prop(Object) step;
  @Prop(Function) deleteStep;

  readonly required = required;

  readonly metrics = [
    "euclidean",
    "manhattan",
    "chebyshev",
    "minkowski",
    "canberra",
    "braycurtis",
    "haversine",
    "mahalanobis",
    "wminkowski",
    "seuclidean",
    "cosine",
    "correlation",
  ];
}
</script>
