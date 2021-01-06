<template>
  <v-container fluid>
    <v-toolbar dense>
      <v-toolbar-title>Edit Model</v-toolbar-title>
      <v-spacer />
      <v-toolbar-items>
        <v-btn @click="cancel" text color="primary">Cancel</v-btn>
        <v-btn @click="reset" text color="primary">Reset</v-btn>
        <v-btn @click="submit" text :disabled="!valid" color="primary">Save</v-btn>
      </v-toolbar-items>
    </v-toolbar>
    <v-card class="mt-4 px-4">
      <v-card-text>
        <v-form v-model="valid" ref="form" lazy-validation>
          <v-text-field label="Name" v-model="name" :rules="nameRules" />
          <v-text-field label="Description" v-model="description" />
        </v-form>
      </v-card-text>
    </v-card>
  </v-container>
</template>

<script lang="ts">
import { Component, Vue } from "vue-property-decorator";
import { required } from "@/utils/validators";
import { modelsModule } from "@/modules/models";
import { IModelUpdate } from "@/modules/models/models";

@Component
export default class EditModel extends Vue {
  readonly modelsContext = modelsModule.context(this.$store);

  readonly nameRules = [required];

  valid = true;
  name = "";
  description: string | null = null;

  get model() {
    return this.modelsContext.getters.getModel(+this.$router.currentRoute.params.id);
  }

  async mounted() {
    await this.modelsContext.actions.getModel(+this.$router.currentRoute.params.id);
    this.reset();
  }

  reset() {
    if (this.$refs.form) {
      (this.$refs.form as any).resetValidation();
    }
    if (this.model) {
      this.name = this.model.name;
      this.description = this.model.description;
    }
  }

  cancel() {
    this.$router.back();
  }

  async submit() {
    if ((this.$refs.form as any).validate() && this.model) {
      const data: IModelUpdate = {
        name: this.name,
        description: this.description,
      };
      await this.modelsContext.actions.updateModel({
        modelId: this.model.id,
        data: data,
      });
      this.$router.back();
    }
  }
}
</script>
