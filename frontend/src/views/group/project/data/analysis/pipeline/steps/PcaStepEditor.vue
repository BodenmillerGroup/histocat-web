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
        PCA
      </v-expansion-panel-header>
      <v-expansion-panel-content>
        <v-card flat>
          <v-radio-group label="SVD solver" v-model="step.svdSolver">
            <v-radio label="ARPACK" value="arpack" />
            <v-radio label="randomized" value="randomized" />
            <v-radio label="auto" value="auto" />
            <v-radio label="LOBPCG" value="lobpcg" />
          </v-radio-group>
        </v-card>
      </v-expansion-panel-content>
    </v-expansion-panel>
  </v-expansion-panels>
</template>

<script lang="ts">
import { Component, Prop, Vue } from "vue-property-decorator";

@Component
export default class PcaStepEditor extends Vue {
  @Prop(Object) step;
  @Prop(Function) deleteStep;
}
</script>
