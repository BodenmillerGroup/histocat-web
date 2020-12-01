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
        UMAP
      </v-expansion-panel-header>
      <v-expansion-panel-content>
        <v-card flat>
          <v-row dense no-gutters>
            <v-col>
              <v-text-field
                label="Minimum distance"
                hint="The effective minimum distance between embedded points. Smaller values will result in a more clustered/clumped embedding where nearby points on the manifold are drawn closer together, while larger values will result on a more even dispersal of points. The value should be set relative to the ``spread`` value, which determines the scale at which embedded points will be spread out. The default of in the `umap-learn` package is 0.1."
                v-model.number="step.minDist"
                type="number"
                min="0"
                :rules="[required]"
                class="text-field"
              />
            </v-col>
            <v-col>
              <v-text-field
                type="number"
                min="0"
                step="0.1"
                hint="The effective scale of embedded points. In combination with `min_dist` this determines how clustered/clumped the embedded points are."
                label="Spread"
                v-model.number="step.spread"
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
  @Prop({ type: Object, required: true }) step;
  @Prop({ type: Function, required: true }) deleteStep;
  @Prop({ type: Function, required: true }) moveUpStep;
  @Prop({ type: Function, required: true }) moveDownStep;

  readonly required = required;
}
</script>

<style scoped>
.text-field {
  margin-right: 20px;
}
</style>
