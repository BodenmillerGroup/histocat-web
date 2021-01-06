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
        </v-form>
      </v-card-text>
    </v-card>
    <v-card class="mt-4 px-4">
      <v-card-title primary-title>
        <div class="text-h5">Model File</div>
      </v-card-title>
      <v-card-text>
        <v-form>
          <v-file-input v-model="file" label="File upload" show-size />
        </v-form>
      </v-card-text>
    </v-card>
  </v-container>
</template>

<script lang="ts">
import { Component, Vue } from "vue-property-decorator";
import { groupModule } from "@/modules/group";
import { required } from "@/utils/validators";
import { IModelCreate } from "@/modules/models/models";
import { modelsModule } from "@/modules/models";

@Component
export default class CreateModel extends Vue {
  readonly groupContext = groupModule.context(this.$store);
  readonly modelsContext = modelsModule.context(this.$store);

  readonly nameRules = [required];

  valid = true;
  name = "";
  description: string | null = null;
  file: File | null = null;

  get activeGroupId() {
    return this.groupContext.getters.activeGroupId;
  }

  reset() {
    this.name = "";
    this.description = null;
    this.file = null;
    (this.$refs.form as any).resetValidation();
  }

  cancel() {
    this.$router.back();
  }

  async submit() {
    if ((this.$refs.form as any).validate() && this.activeGroupId) {
      const data: IModelCreate = {
        name: this.name,
        description: this.description,
      };
      const model = await this.modelsContext.actions.createModel(data);

      if (model && model.id && this.file) {
        const formData = new FormData();
        formData.append("groupId", this.activeGroupId.toString());
        formData.append("file", this.file);
        await this.modelsContext.actions.uploadModelFile({
          modelId: model.id,
          formData: formData,
        });
      }

      this.$router.back();
    }
  }
}
</script>
