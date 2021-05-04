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
        t-SNE
      </v-expansion-panel-header>
      <v-expansion-panel-content>
        <v-card flat>
          <v-row dense no-gutters>
            <v-col>
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
            </v-col>
            <v-col>
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
            </v-col>
            <v-col>
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
export default class TsneStepEditor extends Vue {
  @Prop({ type: Object, required: true }) readonly step;
  @Prop({ type: Function, required: true }) readonly deleteStep;
  @Prop({ type: Function, required: true }) readonly moveUpStep;
  @Prop({ type: Function, required: true }) readonly moveDownStep;

  readonly required = required;
}
</script>

<style scoped>
.text-field {
  margin-right: 20px;
}
</style>
