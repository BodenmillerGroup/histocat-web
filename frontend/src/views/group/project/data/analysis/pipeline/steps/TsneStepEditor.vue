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
        t-SNE
      </v-expansion-panel-header>
      <v-expansion-panel-content>
        <v-card flat>
          <v-row dense no-gutters>
            <v-col cols="4">
              <v-text-field
                label="Perplexity"
                hint="The perplexity is related to the number of nearest neighbors that is used in other manifold learning algorithms. Larger datasets usually require a larger perplexity. Consider selecting a value between 5 and 50. The choice is not extremely critical since t-SNE is quite insensitive to this parameter."
                v-model.number="step.perplexity"
                type="number"
                min="2"
                max="100"
                step="1"
                :rules="[required]"
                class="text-field"
              />
              <v-text-field
                label="Early exaggeration"
                hint="Controls how tight natural clusters in the original space are in the embedded space and how much space will be between them. For larger values, the space between natural clusters will be larger in the embedded space. Again, the choice of this parameter is not very critical. If the cost function increases during initial optimization, the early exaggeration factor or the learning rate might be too high."
                v-model.number="step.earlyExaggeration"
                type="number"
                min="0"
                step="1"
                :rules="[required]"
                class="text-field"
              />
              <v-text-field
                label="Learning rate"
                hint="Note that the R-package “Rtsne” uses a default of 200. The learning rate can be a critical parameter. It should be between 100 and 1000. If the cost function increases during initial optimization, the early exaggeration factor or the learning rate might be too high. If the cost function gets stuck in a bad local minimum increasing the learning rate helps sometimes."
                v-model.number="step.learningRate"
                type="number"
                min="100"
                max="1000"
                step="1"
                :rules="[required]"
                class="text-field"
              />
            </v-col>
            <v-col>
              <VariablesSelector :step="step" />
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
import VariablesSelector from "@/views/group/project/data/analysis/pipeline/steps/VariablesSelector.vue";

@Component({
  components: { VariablesSelector },
})
export default class TsneStepEditor extends Vue {
  @Prop(Object) step;
  @Prop(Function) deleteStep;

  readonly required = required;
}
</script>

<style scoped>
.text-field {
  margin-right: 20px;
}
</style>
