<template>
  <v-container fluid>
    <v-toolbar dense>
      <v-toolbar-title>Add Model</v-toolbar-title>
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
          <v-select :items="applications" v-model="application" label="Application" dense :rules="applicationRules" />
          <v-file-input v-model="file" label="Model file" show-size accept=".zip" :rules="fileRules" />
        </v-form>
      </v-card-text>
    </v-card>
  </v-container>
</template>

<script lang="ts">
import { Component, Vue } from "vue-property-decorator";
import { required } from "@/utils/validators";
import { modelsModule } from "@/modules/models";

@Component
export default class CreateModel extends Vue {
  readonly modelsContext = modelsModule.context(this.$store);

  readonly nameRules = [required];
  readonly fileRules = [required];
  readonly applicationRules = [required];

  readonly applications = ["mesmer", "nuclear", "cytoplasm"];

  valid = true;
  name = "";
  application = "mesmer";
  description = "";
  file: File | null = null;

  reset() {
    this.name = "";
    this.description = "";
    this.application = "mesmer";
    this.file = null;
    (this.$refs.form as any).resetValidation();
  }

  cancel() {
    this.$router.back();
  }

  async submit() {
    if ((this.$refs.form as any).validate() && this.file) {
      const formData = new FormData();
      formData.append("name", this.name);
      formData.append("application", this.application);
      formData.append("description", this.description);
      formData.append("file", this.file);
      await this.modelsContext.actions.createModel(formData);

      this.$router.back();
    }
  }
}
</script>
