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
        Transformation
      </v-expansion-panel-header>
      <v-expansion-panel-content>
        <v-card flat>
          <v-row dense no-gutters>
            <v-col>
              <v-radio-group label="Mode" v-model="step.mode">
                <v-radio label="arcsinh" value="arcsinh" />
                <v-radio label="log1p" value="log1p" />
              </v-radio-group>
            </v-col>
            <v-col>
              <v-text-field
                label="Cofactor"
                :disabled="step.mode !== 'arcsinh'"
                v-model.number="step.cofactor"
                type="number"
                min="0"
                step="1"
                :rules="cofactorRules"
                class="mt-2 text-field"
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
import { positiveNumber, required } from "@/utils/validators";

@Component
export default class TransformationStepEditor extends Vue {
  @Prop({ type: Object, required: true }) step;
  @Prop({ type: Function, required: true }) deleteStep;
  @Prop({ type: Function, required: true }) moveUpStep;
  @Prop({ type: Function, required: true }) moveDownStep;

  readonly cofactorRules = [required, positiveNumber];
}
</script>

<style scoped>
.text-field {
  margin-right: 20px;
}
</style>
