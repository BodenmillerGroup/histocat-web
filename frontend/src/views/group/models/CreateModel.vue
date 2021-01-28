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
          <v-file-input v-model="file" label="Model file" show-size accept=".zip" :rules="fileRules" />
        </v-form>
      </v-card-text>
    </v-card>
  </v-container>
</template>

<script lang="ts">
import { Component, Vue } from "vue-property-decorator";
import { groupModule } from "@/modules/group";
import { required } from "@/utils/validators";
import { modelsModule } from "@/modules/models";

@Component
export default class CreateModel extends Vue {
  readonly groupContext = groupModule.context(this.$store);
  readonly modelsContext = modelsModule.context(this.$store);

  readonly nameRules = [required];
  readonly fileRules = [required];

  valid = true;
  name = "";
  description = "";
  file: File | null = null;

  get activeGroupId() {
    return this.groupContext.getters.activeGroupId;
  }

  reset() {
    this.name = "";
    this.description = "";
    this.file = null;
    (this.$refs.form as any).resetValidation();
  }

  cancel() {
    this.$router.back();
  }

  async submit() {
    if ((this.$refs.form as any).validate() && this.activeGroupId && this.file) {
      const formData = new FormData();
      formData.append("name", this.name);
      formData.append("description", this.description);
      formData.append("file", this.file);
      await this.modelsContext.actions.createModel(formData);

      this.$router.back();
    }
  }
}
</script>
