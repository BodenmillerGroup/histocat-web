<template>
  <v-container fluid>
    <v-toolbar dense>
      <v-toolbar-title>
        Edit Group
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
import {IGroupUpdate} from "@/modules/group/models";

@Component
export default class EditGroup extends Vue {
  readonly groupContext = groupModule.context(this.$store);

  readonly nameRules = [required];

  valid = true;
  name = "";
  description: string | null = null;
  url: string | null = null;
  isOpen = false;

  get group() {
    return this.groupContext.getters.getGroup(+this.$router.currentRoute.params.groupId);
  }

  async mounted() {
    await this.groupContext.actions.getGroup(+this.$router.currentRoute.params.groupId);
    this.reset();
  }

  reset() {
    this.name = "";
    this.description = "";
    this.url = "";
    this.isOpen = false;
    if (this.$refs.form) {
      (this.$refs.form as any).resetValidation();
    }
    if (this.group) {
      this.name = this.group.name;
      this.description = this.group.description;
      this.url = this.group.url;
      this.isOpen = this.group.is_open;
    }
  }

  cancel() {
    this.$router.back();
  }

  async submit() {
    if ((this.$refs.form as any).validate()) {
      const data: IGroupUpdate = {
        name: this.name,
        description: this.description,
        url: this.url,
        is_open: this.isOpen,
      };
      await this.groupContext.actions.updateGroup({
        id: this.group!.id,
        data: data,
      });
      this.$router.back();
    }
  }
}
</script>
