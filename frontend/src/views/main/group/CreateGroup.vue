<template>
  <v-container fluid>
    <v-toolbar dense>
      <v-toolbar-title>
        Create Group
      </v-toolbar-title>
      <v-spacer />
      <v-toolbar-items>
        <v-btn @click="cancel" text color="primary">Cancel</v-btn>
        <v-btn @click="reset" text color="primary">Reset</v-btn>
        <v-btn @click="submit" text :disabled="!valid" color="primary">Save</v-btn>
      </v-toolbar-items>
    </v-toolbar>
    <v-card class="mt-4 px-4">
      <v-card-text>
        <template>
          <v-form v-model="valid" ref="form" lazy-validation>
            <v-text-field label="Name" v-model="name" :rules="nameRules" />
            <v-text-field label="Description" v-model="description" />
            <v-text-field label="URL" v-model="url" />
            <v-checkbox label="Open" v-model="isOpen" />
          </v-form>
        </template>
      </v-card-text>
    </v-card>
  </v-container>
</template>

<script lang="ts">
import { required } from "@/utils/validators";
import { Component, Vue } from "vue-property-decorator";
import { groupModule } from "@/modules/group";
import { IGroupCreate } from "@/modules/group/models";

@Component
export default class CreateGroup extends Vue {
  readonly groupContext = groupModule.context(this.$store);

  readonly nameRules = [required];

  valid = true;
  name = "";
  description: string | null = null;
  url: string | null = null;
  isOpen = false;

  reset() {
    this.name = "";
    this.description = null;
    this.url = null;
    this.isOpen = false;
    (this.$refs.form as any).resetValidation();
  }

  cancel() {
    this.$router.back();
  }

  async submit() {
    if ((this.$refs.form as any).validate()) {
      const params: IGroupCreate = {
        name: this.name,
        description: this.description,
        url: this.url,
        is_open: this.isOpen,
      };
      await this.groupContext.actions.createGroup(params);
      this.$router.back();
    }
  }
}
</script>
