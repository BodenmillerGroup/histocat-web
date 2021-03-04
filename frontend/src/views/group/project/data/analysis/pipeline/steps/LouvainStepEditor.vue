<template>
  <v-expansion-panels :value="0">
    <v-expansion-panel>
      <v-expansion-panel-header disable-icon-rotate>
        <template v-slot:actions>
          <v-tooltip bottom>
            <template v-slot:activator="{ on }">
              <v-btn icon v-on="on" @click.stop="moveUpStep(step)">
                <v-icon color="primary">mdi-arrow-up-bold-outline</v-icon>
              </v-btn>
            </template>
            <span>Move step up</span>
          </v-tooltip>
          <v-tooltip bottom>
            <template v-slot:activator="{ on }">
              <v-btn icon v-on="on" @click.stop="moveDownStep(step)">
                <v-icon color="primary">mdi-arrow-down-bold-outline</v-icon>
              </v-btn>
            </template>
            <span>Move step down</span>
          </v-tooltip>
          <v-tooltip bottom>
            <template v-slot:activator="{ on }">
              <v-btn icon v-on="on" @click.stop="deleteStep(step)">
                <v-icon color="secondary">mdi-delete-outline</v-icon>
              </v-btn>
            </template>
            <span>Delete step</span>
          </v-tooltip>
        </template>
        Louvain clustering
      </v-expansion-panel-header>
      <v-expansion-panel-content>
        <v-card flat>
          <v-row dense>
            <v-col>
              <v-text-field
                label="Resolution"
                hint="A parameter value controlling the coarseness of the clustering. Higher values lead to more clusters."
                v-model.number="step.resolution"
                type="number"
                min="0"
                step="0.1"
                class="mr-4"
                :rules="[required]"
              />
              <v-select dense label="Flavor" :items="flavors" v-model="step.flavor" class="mr-4" />
            </v-col>
            <v-col>
              <v-checkbox
                label="Directed"
                hint="Whether to treat the graph as directed or undirected."
                persistent-hint
                v-model="step.directed"
              />
              <v-checkbox
                label="Use weights"
                hint="Use weights from knn graph."
                persistent-hint
                v-model="step.useWeights"
              />
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
export default class LouvainStepEditor extends Vue {
  @Prop({ type: Object, required: true }) step;
  @Prop({ type: Function, required: true }) deleteStep;
  @Prop({ type: Function, required: true }) moveUpStep;
  @Prop({ type: Function, required: true }) moveDownStep;

  readonly required = required;

  readonly flavors = ["vtraag", "igraph", "rapids"];
}
</script>
