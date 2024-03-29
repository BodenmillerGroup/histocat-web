<template>
  <v-container fluid>
    <v-toolbar dense>
      <v-toolbar-title>Edit Group Member</v-toolbar-title>
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
          <div class="text-subtitle-1">Role</div>
          <v-btn-toggle v-model="role">
            <v-btn small value="100">Admin</v-btn>
            <v-btn small value="10">Standard</v-btn>
            <v-btn small value="0">Guest</v-btn>
          </v-btn-toggle>
          <v-checkbox label="Active" v-model="isActive" hint="Access is permited" />
        </v-form>
      </v-card-text>
    </v-card>
  </v-container>
</template>

<script lang="ts">
import { Component, Vue } from "vue-property-decorator";
import { memberModule } from "@/modules/member";
import { IMemberUpdate } from "@/modules/member/models";

@Component
export default class EditMember extends Vue {
  readonly memberContext = memberModule.context(this.$store);

  valid = true;
  role = "0";
  isActive = false;

  get member() {
    return this.memberContext.getters.getMember(+this.$router.currentRoute.params.id);
  }

  async mounted() {
    await this.memberContext.actions.getMember(+this.$router.currentRoute.params.id);
    this.reset();
  }

  reset() {
    if (this.$refs.form) {
      (this.$refs.form as any).resetValidation();
    }
    if (this.member) {
      this.role = this.member.role.toString();
      this.isActive = this.member.is_active;
    }
  }

  cancel() {
    this.$router.back();
  }

  async submit() {
    if ((this.$refs.form as any).validate() && this.member) {
      const data: IMemberUpdate = {
        role: !this.role || this.role === "" ? 0 : Number(this.role),
        is_active: this.isActive,
      };
      await this.memberContext.actions.updateMember({
        id: this.member.id,
        data: data,
      });
      this.$router.back();
    }
  }
}
</script>
